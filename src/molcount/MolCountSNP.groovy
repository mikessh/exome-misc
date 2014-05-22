package molcount

@Grab(group = 'org.utgenome.thirdparty', module = 'picard', version = '1.86.0')

import net.sf.samtools.AlignmentBlock
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

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
    println "groovy MolCountSNP input(sam or bam) varinats(VCF format) output [read length, default = 100]"
    println "=DESCR"
    println "Estimates the number of initial molecules using the number of read (pairs) with distinct mapping offsets"
    println "Uses these estimates to calculate 'molecular coverage' of each SNP in VCF data"
    System.exit(0)
}

def scriptName = getClass().canonicalName
def inputSamFileName = args[0], inputMutectFileName = args[1], outputFileName = args[2]
int readLength = args.length > 3 ? Integer.parseInt(args[3]) : 200 // offset for overlapping read search by hash

// Build variant database
class VariantMolData {
    String var, ref
    int pos
    def molsWithVar = new HashSet<Integer>(),
        molsWithRef = new HashSet<Integer>(),
        molsWithOther = new HashSet<Integer>()

    int readsWithVar, readsWithRef, readsWithOther

    String getStatsString() {
        molsWithVar.size() + "\t" + molsWithRef.size() + "\t" + molsWithOther.size() + "\t" +
                readsWithVar + "\t" + readsWithRef + "\t" + readsWithOther
    }

    static String getStatsHeader() {
        "MOLECULES_WITH_VARIANT\tMOLECULES_WITH_REF\tMOLECULES_WITH_OTHER\t" +
                "READS_WITH_VARIANT\tREADS_WITH_REF\tREADS_WITH_OTHER"
    }
}

def FIELD_INDICES = [0, 1, 3, 4] // contig, pos, ref, var

def contigPositionVariantMap = new HashMap<String, Map<Integer, List<VariantMolData>>>()
def variantMap = new HashMap<String, VariantMolData>()
println "[${new Date()} $scriptName] Loading variants from $inputMutectFileName.."
new File(inputMutectFileName).withReader { reader ->
    while (reader.readLine().startsWith("#")) {
        //
    }

    String line
    while ((line = reader.readLine()) != null) {
        def splitLine = line.split("\t")
        def contig = splitLine[FIELD_INDICES[0]]
        int pos = Integer.parseInt(splitLine[FIELD_INDICES[1]])
        def ref = splitLine[FIELD_INDICES[2]]
        def var = splitLine[FIELD_INDICES[3]]

        def variantMolData = new VariantMolData(var: var, ref: ref, pos: pos)

        def positionVariantMap = contigPositionVariantMap.get(contig)

        if (positionVariantMap == null)
            contigPositionVariantMap.put(contig, positionVariantMap = new HashMap<Integer, List<VariantMolData>>())

        for (int i = -readLength; i <= readLength; i++) {
            def variantMolCountDataList = positionVariantMap.get(pos + i)
            if (variantMolCountDataList == null)
                positionVariantMap.put(pos + i, variantMolCountDataList = new ArrayList<VariantMolData>())
            variantMolCountDataList.add(variantMolData)
        }

        variantMap.put(contig + ":" + pos + ">" + var, variantMolData)
    }
}
println "[${new Date()} $scriptName] ${variantMap.size()} variants loaded"

def getBaseByRefPos = { int pos, SAMRecord record ->
    // Find a position in read by reference coordinate
    // We use contigously aligned block coordinates

    for (AlignmentBlock alignmentBlock : record.alignmentBlocks) {
        int delta = pos - alignmentBlock.referenceStart
        if (delta > 0 && delta < alignmentBlock.length - 1) {
            return record.readString[alignmentBlock.readStart + delta - 1]
        }
    }
    'N'
}

// Go through SAM file
def samReader = new SAMFileReader(new File(inputSamFileName))
int nReads = 0, nReadsVar = 0, nReadsOther = 0, nReadsRef = 0, nReadsUndef = 0
def unusedContigs = new HashSet<String>()
println "[${new Date()} $scriptName] Going through $inputSamFileName.."
samReader.each { SAMRecord record ->
    if (!record.readUnmappedFlag) { // Mapped read
        def positionVariantMap = contigPositionVariantMap.get(record.referenceName)

        if (positionVariantMap != null) { // Chromosome of interest
            def variantMolCountDataList = positionVariantMap.get(record.alignmentStart)
            if (variantMolCountDataList != null) { // Any SNPs in vicinity?
                variantMolCountDataList.each { variantMolCountData ->
                    if (variantMolCountData.pos >= record.alignmentStart &&
                            variantMolCountData.pos <= record.alignmentEnd) {
                        String readVariant = getBaseByRefPos(variantMolCountData.pos, record)
                        if (readVariant != 'N') {
                            if (readVariant == variantMolCountData.ref) {
                                variantMolCountData.molsWithRef.add(
                                        (record.alignmentStart +
                                                (record.mateUnmappedFlag ? 0 : record.mateAlignmentStart) +
                                                readLength) / 2 as int)
                                variantMolCountData.readsWithRef++
                                nReadsRef++
                            } else if (readVariant == variantMolCountData.var) {
                                variantMolCountData.molsWithVar.add(
                                        (record.alignmentStart +
                                                (record.mateUnmappedFlag ? 0 : record.mateAlignmentStart) +
                                                readLength) / 2 as int)
                                variantMolCountData.readsWithVar++
                                nReadsVar++
                            } else {
                                variantMolCountData.molsWithOther.add(
                                        (record.alignmentStart +
                                                (record.mateUnmappedFlag ? 0 : record.mateAlignmentStart) +
                                                readLength) / 2 as int)
                                variantMolCountData.readsWithOther++
                                nReadsOther++
                            }
                        } else
                            nReadsUndef++
                    }
                }
            }
        }

        if (++nReads % 1000000 == 0)
            println "[${new Date()} $scriptName] $nReads reads processed, currently " +
                    "at ${record.referenceName}:${record.alignmentStart}. " +
                    "READS: VAR_ALIGN=$nReadsVar REF_ALIGN=$nReadsRef OTHER_ALIGN=$nReadsOther " +
                    "UNDEF_ALIGN=$nReadsUndef"
    }
}
println "[${new Date()} $scriptName] Finished, $nReads reads processed"
println "[${new Date()} $scriptName] Unused contigs: $unusedContigs"

println "[${new Date()} $scriptName] Writing outut to $outputFileName.."
new File(outputFileName).withPrintWriter { pw ->
    new File(inputMutectFileName).withReader { reader ->
        String line
        while ((line = reader.readLine()).startsWith("#")) {
            pw.println(line)
        }
        pw.println(line + "\t" + VariantMolData.statsHeader)

        while ((line = reader.readLine()) != null) {
            def splitLine = line.split("\t")
            def contig = splitLine[FIELD_INDICES[0]]
            int pos = Integer.parseInt(splitLine[FIELD_INDICES[1]])
            def var = splitLine[FIELD_INDICES[3]]

            def variantMolData = variantMap.get(contig + ":" + pos + ">" + var)

            pw.println(line + "\t" + variantMolData.statsString)
        }
    }
}