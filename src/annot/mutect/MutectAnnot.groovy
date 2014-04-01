package annot.mutect

/**
 Copyright 2013 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

// CLI
def cli = new CliBuilder(usage: 'MutectAnnot.groovy [options] mutect_output1[,mutect_output2,...] [path to folder with refGene, COSMIC and dbSNP data if not in script dir]')
cli.h('Display help')
cli.i(longOpt: 'intersect', 'Report intersections between vcf files (adds columns with tumor allele frequencies in all datasets)')
cli.u('Append URL links to output')
cli.e('Extended reporting: include all transcript isoforms. Default: only most representative isoform is used, with segment hierarchy Exon>UTR>Intron then selected by CDS length')
def opt = cli.parse(args)
if (opt.arguments().size() < 1 || opt.h) {
    cli.usage()
    System.exit(1)
}

def BIN_SZ = 100000
def header = ""
boolean simple = !opt.e, outputUrl = opt.u, intersect = opt.i
def vcfFiles = opt.arguments()[0].split(",").collect { new File(it) }
def dataPath = opt.arguments().size() < 2 ? new File(getClass().protectionDomain.codeSource.location.path).parent.replaceAll("%20", " ") : opt.arguments()[1]

class Alphabet {
    static char compl(char c) {
        switch (c) {
            case 'A': return 'T'
            case 'T': return 'A'
            case 'G': return 'C'
            case 'C': return 'G'
            case 'a': return 't'
            case 't': return 'a'
            case 'g': return 'c'
            case 'c': return 'g'
            default: return 'N'
        }
    }

    static String codon2aa(String codon) {
        switch (codon.toUpperCase()) {
            case 'TTT': return 'Phe'
            case 'TTC': return 'Phe'
            case 'TTA': return 'Leu'
            case 'TTG': return 'Leu'
            case 'CTT': return 'Leu'
            case 'CTC': return 'Leu'
            case 'CTA': return 'Leu'
            case 'CTG': return 'Leu'
            case 'ATT': return 'Ile'
            case 'ATC': return 'Ile'
            case 'ATA': return 'Ile'
            case 'ATG': return 'Met'
            case 'GTT': return 'Val'
            case 'GTC': return 'Val'
            case 'GTA': return 'Val'
            case 'GTG': return 'Val'
            case 'TCT': return 'Ser'
            case 'TCC': return 'Ser'
            case 'TCA': return 'Ser'
            case 'TCG': return 'Ser'
            case 'CCT': return 'Pro'
            case 'CCC': return 'Pro'
            case 'CCA': return 'Pro'
            case 'CCG': return 'Pro'
            case 'ACT': return 'Thr'
            case 'ACC': return 'Thr'
            case 'ACA': return 'Thr'
            case 'ACG': return 'Thr'
            case 'GCT': return 'Ala'
            case 'GCC': return 'Ala'
            case 'GCA': return 'Ala'
            case 'GCG': return 'Ala'
            case 'TAT': return 'Tyr'
            case 'TAC': return 'Tyr'
            case 'TAA': return '***'
            case 'TAG': return '***'
            case 'CAT': return 'His'
            case 'CAC': return 'His'
            case 'CAA': return 'Gln'
            case 'CAG': return 'Gln'
            case 'AAT': return 'Asn'
            case 'AAC': return 'Asn'
            case 'AAA': return 'Lys'
            case 'AAG': return 'Lys'
            case 'GAT': return 'Asp'
            case 'GAC': return 'Asp'
            case 'GAA': return 'Glu'
            case 'GAG': return 'Glu'
            case 'TGT': return 'Cys'
            case 'TGC': return 'Cys'
            case 'TGA': return '***'
            case 'TGG': return 'Trp'
            case 'CGT': return 'Arg'
            case 'CGC': return 'Arg'
            case 'CGA': return 'Arg'
            case 'CGG': return 'Arg'
            case 'AGT': return 'Ser'
            case 'AGC': return 'Ser'
            case 'AGA': return 'Arg'
            case 'AGG': return 'Arg'
            case 'GGT': return 'Gly'
            case 'GGC': return 'Gly'
            case 'GGA': return 'Gly'
            case 'GGG': return 'Gly'
            default: return '???'
        }
    }

    static String translate(String ntSeq) {
        return ntSeq.split(" ").collect { codon2aa(it) }.join("")
    }

    static String revCompl(String s) {
        def sb = new StringBuilder()
        for (int i = s.length() - 1; i >= 0; i--)
            sb.append(compl(s.charAt(i)))
        sb.toString()
    }
}

// Pre-load refGene
println "Pre-loading RefGene data"
class GeneEntry {
    String geneName, transcriptName
    String chr
    int start, end, cdsStart, cdsEnd
    def exons = new ArrayList<ExonEntry>()
    def exonStarts = new ArrayList<Integer>(), exonEnds = new ArrayList<Integer>()
    def codingSequence = ""
    int strand

    void splice() {
        exons.sort { it.id }
        codingSequence = exons.collect { it.sequence }.join("")
    }

    def getCodons(int coord, int exonId, VcfEntry vcfEntry) {
        def tail = 0
        for (int i = 0; i < exonId; i++)
            tail += exons[i].length()

        tail += strand > 0 ? coord - exons[exonId].start : exons[exonId].end - coord
        def frame = tail % 3, codon = tail - frame

        def ref, var
        ref = codingSequence.substring(codon, codon + 3)
        var = new StringBuilder(ref)
        var.setCharAt(frame, strand > 0 ? (char) vcfEntry.var : Alphabet.compl((char) vcfEntry.var))

        var = var.toString()
        String aaref = Alphabet.translate(ref), aavar = Alphabet.translate(var)
        return [ref, var, aaref + (codon / 3 + 1).toString() + aavar,
                vcfEntry.var.length() > 1 ? (((vcfEntry.var.length() - vcfEntry.ref.length()) % 3 == 0) ? "Frameshift" : ".") :
                        ((aaref != aavar) ? "Missense" : ".")]
    }

    def getSegment(VcfEntry vcfEntry) {
        def coord = vcfEntry.coord

        def eid = -1
        for (int i = 0; i < exons.size(); i++)
            if (coord >= exons[i].start && coord < exons[i].end)
                eid = i

        for (int i = 0; i < exonStarts.size(); i++) { // This is used to get correct exon number, as utr exons are missing in .fa
            if (coord >= exonStarts[i] && coord <= exonEnds[i]) {
                if (coord < cdsStart || coord > cdsEnd)
                    return coord > cdsEnd ? ["3'UTR", -1, coord - cdsEnd] : ["5'UTR", -1, coord - cdsStart]
                else
                    return ["Exon${strand > 0 ? i + 1 : exonStarts.size() - i}",
                            eid, strand > 0 ? coord - exonStarts[i] : exonEnds[i] - coord]
            } else if (coord < exonStarts[i]) {
                def mid = (exonEnds[i - 1] - exonStarts[i]) / 2
                coord = (coord - exonEnds[i - 1]) > mid ? coord - exonStarts[i] : coord - exonEnds[i]
                return ["Intron${strand > 0 ? i : exonStarts.size() - i}", eid, strand * coord]
            }
        }
    }
}

class ExonEntry {
    int id
    int start, end
    String sequence

    int length() {
        end - start + 1
    }
}

def geneMapByBin = new HashMap<String, List<GeneEntry>>(), geneMapByName = new HashMap<String, GeneEntry>()

// Fetch gene coords
// #name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	name2   exonStarts  exonEnds
// 0        1       2       3       4       5           6       7           8           9       10
new File(dataPath + "/refseq_coords.txt").splitEachLine("\t") {
    splitLine ->
        if (!splitLine[0].startsWith("#")) {
            def chrBinStart = splitLine[1] + ":" + (int) ((int) Integer.parseInt(splitLine[3]) / (int) BIN_SZ),
                chrBinEnd = splitLine[1] + ":" + (int) ((int) Integer.parseInt(splitLine[4]) / (int) BIN_SZ)
            def geneEntry = new GeneEntry(transcriptName: splitLine[0],
                    chr: splitLine[1], strand: splitLine[2] == "+" ? 1 : -1, start: Integer.parseInt(splitLine[3]), end: Integer.parseInt(splitLine[4]),
                    cdsStart: Integer.parseInt(splitLine[5]), cdsEnd: Integer.parseInt(splitLine[6]),
                    exonStarts: splitLine[7].split(",").collect { Integer.parseInt(it) }.asList(),
                    exonEnds: splitLine[8].split(",").collect { Integer.parseInt(it) }.asList(),
                    geneName: splitLine[9]
            )
            geneMapByName.put(splitLine[0], geneEntry)
            [chrBinStart, chrBinEnd].each { bin ->
                def geneList = geneMapByBin.get(bin)
                if (geneList == null)
                    geneList = new ArrayList<GeneEntry>()
                geneList.add(geneEntry)
                geneMapByBin.put(bin, geneList)
            }
        }
}

// Fetch sequences
def reader = new BufferedReader(new FileReader(new File(dataPath + "/refseq.fa")))
def line, seq = ""
def prevHeaderData
while ((line = reader.readLine()) != null) {
    if (line.startsWith(">")) {
        // >hg19_refGene_NM_032291_0 range=chr1:67000042-67000051 5'pad=0 3'pad=0 strand=+ repeatMasking=none
        def headerData = (line =~ />hg19_refGene_(N._\d+)_(\d+) range=chr.+:(\d+)-(\d+)/)
        headerData = headerData.size() > 0 ? headerData[0] : null
        if (prevHeaderData != null) {
            def geneData = geneMapByName.get(prevHeaderData[1])
            if (geneData != null) {
                geneData.exons.add(new ExonEntry(id: Integer.parseInt(prevHeaderData[2]),
                        start: Integer.parseInt(prevHeaderData[3]), end: Integer.parseInt(prevHeaderData[4]),
                        sequence: seq))
            }
        }
        prevHeaderData = headerData
        seq = ""
    } else {
        seq += line.replaceAll("[atgc]", "").trim()
    }
}

println "Finalizing gene structure data..."
int coding = 0
geneMapByName.values().each {
    it.splice()
    if (it.codingSequence.length() > 0) coding++
}
println "Of ${geneMapByName.values().size()} transcripts ${coding} are coding and will be used further"

// Read in VCFs
class VcfEntry {
    String data, chrbin, chrcoord
    int coord
    String var, ref
    boolean modified = false
}

println "Pre-reading MUTECT files (for intersection)"
def intersectionMaps = intersect ? vcfFiles.collect { file ->
    def map = new HashMap<String, Double>()
    file.splitEachLine("\t") { splitLine ->
        if (!splitLine[0].startsWith("#") && splitLine[1].isDouble())
            map.put(splitLine[0] + ":" + splitLine[1] + ":" + splitLine[3].toUpperCase() + ">" + splitLine[4].toUpperCase(),
                    Double.parseDouble(splitLine[17]))
    }
    map
} : null

vcfFiles.each { file ->
    def vcfEntries = new HashMap<String, VcfEntry>()
    println "Reading MUTECT entries from $file.name"
    file.splitEachLine("\t") { splitLine ->
        if (!splitLine[0].startsWith("#") && splitLine[1].isDouble()) {
            if (splitLine[-1] != "REJECT") {
                def data = [splitLine[0..4], splitLine[16..17]].flatten().join("\t")

                def id = splitLine[0] + ":" + splitLine[1] + ":" + splitLine[3].toUpperCase() + ">" + splitLine[4].toUpperCase()
                splitLine[4].split(",").each { var ->
                    vcfEntries.put(id, new VcfEntry(data: data,
                            chrbin: splitLine[0] + ":" + (int) ((int) Integer.parseInt(splitLine[1]) / (int) BIN_SZ),
                            chrcoord: splitLine[0] + ":" + splitLine[1],
                            coord: Integer.parseInt(splitLine[1]),
                            ref: splitLine[3],
                            var: var
                    ))
                }
            }
        }
    }
    println "${vcfEntries.size()} entries read"
    header = ["CHROM\tPOS\tCONTEXT\tREF\tALT\tT_LOD_FSTAR\tT_FREQ"].flatten().collect().join("\t")

    def revId = { String id ->
        id.split(":").collect {
            it.contains(">") ? it.split(">").collect { Alphabet.revCompl(it) }.join(">") : it
        }.join(":")
    }

    // COSMIC
    println "Intersecting with COSMIC"
    int n = 0, k = 0
    header += "\tCOSM_ID\tCOSM_URL"
    new File(dataPath + "/cosmic.txt").splitEachLine("\t") { splitLine ->
        def cosm = splitLine[0]
        [splitLine[1], revId(splitLine[1])].each { id ->
            if (vcfEntries.containsKey(id)) {
                def entry = vcfEntries.get(id)   w
                if (!entry.modified) {
                    entry.data += "\t${cosm}\t" + (outputUrl ? "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=${cosm}" : '.')
                    entry.modified = true
                    k++
                }
            }
        }

        n++
        if (n % 100000 == 0)
            println "${n} COSM records scanned, ${k} intersections so far"
    }
    println "Finished ${n} COSM records scanned, ${k} intersections"
    vcfEntries.values().each {
        if (!it.modified)
            it.data += "\t.\t."

        it.modified = false
    }

    // SNP
    println "Intersecting with dbSNP"
    n = 0
    k = 0
    header += "\tSNP_ID\tSNP_URL"
    new File(dataPath + "/dbsnp.txt").splitEachLine("\t") { splitLine ->
        def rs = splitLine[0]
        [splitLine[1], revId(splitLine[1])].each { id ->
            if (vcfEntries.containsKey(id)) {
                def entry = vcfEntries.get(id)
                if (!entry.modified) {
                    entry.data += "\t${rs}\t" + (outputUrl ? "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=${rs}" : '.')
                    entry.modified = true
                    k++
                }
            }
        }
        n++
        if (n % 1000000 == 0)
            println "${n} RS records scanned, ${k} intersections so far"
    }
    println "Finished ${n} RS records scanned, ${k} intersections"
    vcfEntries.values().each {
        if (!it.modified)
            it.data += "\t.\t."

        it.modified = false
    }

    // REFGENE
    println "Intersecting with REFGENE"
    println "Loading gene structure data"
    header += "\tGENE_ID\tGENECARDS\tTRANSCRIPT_ID\tSEGMENT\tCOORD\tREF_CODON\tVAR_CODON\tAA_SEQ\tTYPE\tUCSC_URL"
    vcfEntries.values().each { vcfEntry ->
        def ucscUrl = outputUrl ? "http://genome.ucsc.edu/cgi-bin/hgTracks?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg19&position=${vcfEntry.chrcoord}" : '.'

        def data = ".\t.\t.\t.\t.\t.\t.\t.\t.\t${ucscUrl}"
        if (geneMapByBin.containsKey(vcfEntry.chrbin)) {
            def transcripts = geneMapByBin.get(vcfEntry.chrbin).findAll { // Find all -> select first good var (with a codon)
                it.start <= vcfEntry.coord && it.end >= vcfEntry.coord && it.codingSequence.length() > 0
            }

            if (transcripts.size() > 0) {
                def trSegPairs = transcripts.collect { [it, it.getSegment(vcfEntry)] }
                if (simple) {
                    def bestEntries = trSegPairs.findAll { it[1][0].startsWith("Exon") }
                    if (bestEntries.size() == 0)
                        bestEntries = trSegPairs.findAll { it[1][0].contains("UTR") }
                    if (bestEntries.size() == 0)
                        bestEntries = trSegPairs
                    trSegPairs = [bestEntries.max { it[0].codingSequence.size() }]
                }
                def codons = ['.', '.', '.', '.']
                def aaList = [], types = []
                trSegPairs.each {
                    if (it[1][1] >= 0) {
                        codons = (it[0] as GeneEntry).getCodons(vcfEntry.coord, it[1][1], vcfEntry)
                        aaList.add(codons[2])
                        types.add(codons[3])
                    } else {
                        aaList.add('.')
                        types.add('.')
                    }
                }
                codons[2] = aaList.join(';')
                codons[3] = types.join(';')
                def geneName = (trSegPairs[0][0] as GeneEntry).geneName
                data = [geneName,
                        outputUrl ? "http://www.genecards.org/cgi-bin/carddisp.pl?gene=$geneName" : '.',
                        trSegPairs.collect { (it[0] as GeneEntry).transcriptName }.join(";"),
                        trSegPairs.collect { it[1][0] }.join(";"),
                        trSegPairs.collect { it[1][2] }.join(";"),
                        codons,
                        ucscUrl].flatten().join("\t")
            }
        }
        vcfEntry.data += "\t" + data
    }

    // INTERSECT
    if (intersect) {
        header += "\t" + vcfFiles.collect { it.name + "_TUMOR_F" }.join("\t")
        vcfEntries.each { vcfEntry ->
            vcfEntry.value.data += "\t" + intersectionMaps.collect { it.get(vcfEntry.key) ?: "." }.join("\t")
        }
    }

    // OUT
    println "Finished"
    new File((file.parent ? "$file.parent/" : "") + file.name.replaceFirst(~/\.[^\.]+$/, '.annot.txt')).withPrintWriter {
        pw ->
            pw.println(header)
            vcfEntries.values().each { pw.println(it.data) }
    }
}