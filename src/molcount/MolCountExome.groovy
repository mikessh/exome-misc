package molcount

@Grapes(
        @Grab(group = 'org.utgenome.thirdparty', module = 'picard', version = '1.86.0')
)

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

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

if (args.length < 3) {
    println "=USAGE"
    println "groovy MolCountExome input(sam or bam) exome_ref(refseq format) output [read length, default = 100]"
    println "=DESCR"
    println "Estimates the number of initial molecules using the number of read (pairs) with distinct mapping offsets"
    println "Uses these estimates to calculate 'molecular coverage' of exome data"
    println "=NOTE"
    println "RefSeq format table should have the following columns (by default in UCSC GB table):"
    println "#name chrom strand txStart txEnd cdsStart cdsEnd exonStarts exonEnds name2"
    System.exit(0)
}

def inputFileName = args[0], inputExonMapFileName = args[1], outputFileName = args[2]
def readLength = args.length > 3 ? args[3] : 100

def exon2Data = new HashMap<String, int[]>()
def chr2Bin2Exons = new HashMap<String, Map<Integer, Set<String>>>()
int binSz = 1000000, exonRegionTolerance = 200
println "[${new Date()} MOLCOUNT] Reading exome data from $inputExonMapFileName"
new File(inputExonMapFileName).splitEachLine("\t") { line ->
    def geneName = line[9] + "\t" + line[0]
    def chr = line[1]

    def exonStarts = line[7].split(",").collect { Integer.parseInt(it) },
        exonEnds = line[8].split(",").collect { Integer.parseInt(it) }

    for (int i = 0; i < exonStarts.size(); i++) {
        def exonSignature = geneName + "\t" + (i + 1) + "\t" + chr

        // Record exon regions
        int[] exonData = new int[4] // 1,2 - start & end; 3,4 - molecules & reads (blank)
        exonData[0] = exonStarts[i] - exonRegionTolerance
        exonData[1] = exonEnds[i] + exonRegionTolerance
        exon2Data.put(exonSignature, exonData) // just store exon data by signature

        // Get or create bin set on chr
        def bin2Exons = chr2Bin2Exons.get(chr)
        if (bin2Exons == null)
            chr2Bin2Exons.put(chr, bin2Exons = new HashMap<Integer, Set<String>>())

        // Put exon link to parent bin
        def bins = [(int) (exonData[0] / binSz), (int) (exonData[1] / binSz)]
        bins.each { bin -> // could span two bins (very rare)
            def exonsList = bin2Exons.get(bin)
            if (exonsList == null)
                bin2Exons.put(bin, exonsList = new HashSet<String>())
            exonsList.add(exonSignature)
        }
    }
}

// Read SAM/BAM
def samReader = new SAMFileReader(new File(inputFileName))
println "[${new Date()} MOLCOUNT] Reading SAM file $inputFileName"
def chr2ReadPos2Count = new HashMap<String, Map<Integer, Integer>>()

// Populate map with seq names
samReader.fileHeader.sequenceDictionary.sequences.each {
    chr2ReadPos2Count.put(it.sequenceName, new HashMap<Integer, Integer>())
}

// Count frequencies of read positions
int nReads = 0
samReader.each { SAMRecord record ->
    def readPos2Count = chr2ReadPos2Count.get(record.referenceName)
    if (readPos2Count != null) { // chromosome of interest
        int pos = (record.alignmentStart +
                record.mateUnmappedFlag ? record.alignmentEnd : record.mateAlignmentStart +
                readLength) / 2

        readPos2Count.put(pos, (readPos2Count.get(pos) ?: 0) + 1)
    }
    if (++nReads % 1000000 == 0)
        println "[${new Date()} MOLCOUNT] $nReads reads processed, currently " +
                "at ${record.referenceName}:${record.alignmentStart}"
}

// Summarize for exons
println "[${new Date()} MOLCOUNT] Calculating read/molecule counts for exome"
nReads = 0
int nMoleculesExome = 0, nReadsExome = 0
chr2ReadPos2Count.each { chrEntry ->
    def chr = chrEntry.key
    def bin2Exons = chr2Bin2Exons.get(chr)
    if (bin2Exons != null) {       // have exons at this chr
        chrEntry.value.each { coordEntry ->
            int coord = coordEntry.key
            int count = coordEntry.value

            // fetch exon
            def exonList = bin2Exons.get((int) (coord / binSz))
            if (exonList != null) { // have exons in this bin
                exonList.each { exonSignature ->
                    int[] exonData = exon2Data.get(exonSignature)
                    if (coord >= exonData[0] && coord <= exonData[1]) {
                        exonData[2] = exonData[2] + 1
                        exonData[3] = exonData[3] + count
                        nMoleculesExome++
                        nReadsExome += count
                    }
                }
            }

            if (++nReads % 1000000 == 0)
                println "[${new Date()} MOLCOUNT] $nReads reads processed, currently at ${chr}:${coord}, " +
                        "reads/molecules intersected: $nMoleculesExome/$nReadsExome"
        }
    }
}

println "[${new Date()} MOLCOUNT] Finished, $nReads reads processed, " +
        "reads/molecules intersected: $nMoleculesExome/$nReadsExome"
new File(outputFileName).withPrintWriter { pw ->
    pw.println("#GENE_NAME\tTRANSCRIPT_ID\tEXON_ID\tCHROM\tSTART\tEND\tMOL_COUNT\tREAD_COUNT")
    exon2Data.each {
        pw.println([it.key, it.value].flatten().join("\t"))
    }
}