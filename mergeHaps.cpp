//Input: 2 bam-files sorted by read-name containing same reads lifted over to reference genome after
//being aligned to the individuals two haplotypes.
//Output: 4 bam-files, one containing the primary alignment for each read pair, 2 containing all haplotype specific alignments and one containing unique non-primary alignments.

#include <stdlib.h>
#include <map>
#include <algorithm>
#include <seqan/bam_io.h>

using namespace std;
using namespace seqan;

struct ReadPairInfo {
    CharString readName;
    int pos1;
    int pos2;
    unsigned hapNum;
    unsigned alignScore;
} ;

bool operator<(const ReadPairInfo & Left, const ReadPairInfo & Right)
{
    if (Left.hapNum < Right.hapNum)
        return true;
    if (Left.hapNum > Right.hapNum)
        return false;
    if (Left.pos1 < Right.pos1)
        return true;
    if (Left.pos1 > Right.pos1)
        return false;
    if (Left.pos2 < Right.pos2)
        return true;
    if (Left.pos2 > Right.pos2)
        return false;
    if (Left.alignScore < Right.alignScore)
        return true;
    if (Left.alignScore > Right.alignScore)
        return false;
    if (Left.readName < Right.readName)
        return true;
    return false;
}

bool operator ==(const ReadPairInfo & Left, const ReadPairInfo & Right)
{
    return Left.readName == Right.readName && Left.pos1 == Right.pos1 && Left.pos2 == Right.pos2 && Left.hapNum == Right.hapNum && Left.alignScore == Right.alignScore;
}

void printStatus(const char * message)
{
    // Get the current date and time.
    char timestamp[80];
    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);
    strftime(timestamp, sizeof(timestamp), "[mergeHaps %Y-%m-%d %X] ", &tstruct);

    // Print time and message.
    std::cerr << timestamp << message << std::endl;
}

void printStatus(std::ostringstream & message)
{
    std::string msg = message.str();
    printStatus(toCString(msg));
}

unsigned findAlignmentScore (CharString& tags)
{
    BamTagsDict tagsDict(tags);
    unsigned tagIdx = 0;
    unsigned alignmentScore = 0;
    if (!findTagKey(tagIdx, tagsDict, "AS"))
    {
        cout << "Can't find alignment score tag, twilight zone mode activated.\n";
    }
    else
    {
        if (!extractTagValue(alignmentScore, tagsDict, tagIdx))
            cout << "Couldn't extract alignment score from tagsDict, engage rifles.\n";
    }
    return alignmentScore;
}

void addToMap(map<ReadPairInfo, Pair<BamAlignmentRecord> >& readPairInfoToReads, BamAlignmentRecord& record, unsigned hapNum)
{
    ReadPairInfo readInfo;
    readInfo.readName = record.qName;
    readInfo.hapNum = hapNum;
    readInfo.alignScore = findAlignmentScore(record.tags);
    if (hasFlagFirst(record))
    {
        readInfo.pos1 = record.beginPos;
        readInfo.pos2 = record.pNext;
        readPairInfoToReads[readInfo].i1 = record;
    }
    else
    {
        readInfo.pos1 = record.pNext;
        readInfo.pos2 = record.beginPos;
        readPairInfoToReads[readInfo].i2 = record;
    }
}

unsigned computeCombinedScore(Pair<BamAlignmentRecord>& readPair)
{
    unsigned combinedScore = 0;
    if (!(hasFlagUnmapped(readPair.i1) && hasFlagUnmapped(readPair.i2)))
    {
        if (readPair.i1.qName != "")
            combinedScore += findAlignmentScore(readPair.i1.tags);
        if (readPair.i2.qName != "")
            combinedScore += findAlignmentScore(readPair.i2.tags);
    }
    return combinedScore;
}

void addTagToReads(Pair<BamAlignmentRecord>& winner, unsigned hapNum)
{
    BamTagsDict tagsDict1stRead(winner.i1.tags);
    BamTagsDict tagsDict2ndRead(winner.i2.tags);
    setTagValue(tagsDict1stRead, "XH", hapNum);
    setTagValue(tagsDict2ndRead, "XH", hapNum);
}

bool alignmentsAreEqual(Pair<BamAlignmentRecord>& hap1reads, Pair<BamAlignmentRecord>& hap2reads)
{
    if (hap1reads.i1.beginPos != hap2reads.i1.beginPos)
        return false;
    if (hap1reads.i2.beginPos != hap2reads.i2.beginPos)
        return false;
    if (hap1reads.i1.cigar != hap2reads.i1.cigar)
        return false;
    if (hap1reads.i2.cigar != hap2reads.i2.cigar)
        return false;
    if (computeCombinedScore(hap1reads) != computeCombinedScore(hap2reads))
        return false;
    return true;
}

Pair<bool, String<Pair<BamAlignmentRecord> > > findUniqueAlignments(map<ReadPairInfo, Pair<BamAlignmentRecord> >& readPairInfoToReads, BamFileOut& hap1specific, BamFileOut& hap2specific)
{
    bool allUnaligned = true;
    String<Pair<BamAlignmentRecord> > readsToPrint;
    Pair<bool, String<Pair<BamAlignmentRecord> > > returnValue;
    ReadPairInfo otherHapInfo;
    vector<ReadPairInfo> keysToIgnore;
    map<ReadPairInfo, Pair<BamAlignmentRecord> >::const_iterator itEnd = readPairInfoToReads.end();
    for(map<ReadPairInfo, Pair<BamAlignmentRecord> >::iterator it = readPairInfoToReads.begin(); it != itEnd; ++it)
    {
        //Check if an equivalent alignment with hap=0 has already been added to list of reads to print.
        if (keysToIgnore.size()>0)
        {
              //If it has been we continue without processing the read-pair
              if (find(keysToIgnore.begin(),keysToIgnore.end(),it->first) != keysToIgnore.end())
                    continue;
        }
        //Look for an equivalent alignment coming from the other haplotype.
        otherHapInfo = it->first;
        if (otherHapInfo.hapNum == 1)
              otherHapInfo.hapNum = 2;
        else
              otherHapInfo.hapNum = 1;
        if (readPairInfoToReads.count(otherHapInfo)!= 0)
        {
              //make absolutely sure that alignments are equal
              if (alignmentsAreEqual(it->second, readPairInfoToReads[otherHapInfo]))
              {
                    //If they are, we set haplotype of current readpair to 0 and add key of equivalent alignment to an ignore list.
                    addTagToReads(it->second,0);
                    keysToIgnore.push_back(otherHapInfo);
              }
              else
              {
                  addTagToReads(it->second, it->first.hapNum);
                  //Write to correct haplotype specific bam file if reads are not both unaligned
                  if (!(hasFlagUnmapped(it->second.i1) && hasFlagUnmapped(it->second.i2)))
                  {
                      if (it->second.i1.qName != "")
                      {
                          if (it->first.hapNum==1)
                              writeRecord(hap1specific,it->second.i1);
                          else
                              writeRecord(hap2specific,it->second.i1);
                      }
                      if (it->second.i2.qName != "")
                      {
                          if (it->first.hapNum==1)
                              writeRecord(hap1specific,it->second.i2);
                          else
                              writeRecord(hap2specific,it->second.i2);
                      }
                  }
              }
        }
        else
        {
            addTagToReads(it->second, it->first.hapNum);
            //Write to correct haplotype specific bam file if reads are not both unaligned
            if (!(hasFlagUnmapped(it->second.i1) && hasFlagUnmapped(it->second.i2)))
            {
                if (it->second.i1.qName != "")
                {
                    if (it->first.hapNum==1)
                        writeRecord(hap1specific,it->second.i1);
                    else
                        writeRecord(hap2specific,it->second.i1);
                }
                if (it->second.i2.qName != "")
                {
                    if (it->first.hapNum==1)
                        writeRecord(hap1specific,it->second.i2);
                    else
                        writeRecord(hap2specific,it->second.i2);
                }
            }
        }
        if (!(hasFlagUnmapped(it->second.i1) && hasFlagUnmapped(it->second.i2)))
            allUnaligned = false;
        appendValue(readsToPrint, it->second);
    }
    returnValue.i1 = allUnaligned;
    returnValue.i2 = readsToPrint;
    return returnValue;
}

unsigned setPrimarySecondaryFlags(String<Pair<BamAlignmentRecord> >& readsToPrint, BamFileOut& mergedHaps, BamFileOut& equalSecondary)
{
    unsigned maxIdx=0, maxScore = 0, currScore;
    String<unsigned> evenReads;
    for (unsigned i=0; i<length(readsToPrint); ++i)
    {
        if (!(hasFlagUnmapped(readsToPrint[i].i1) && hasFlagUnmapped(readsToPrint[i].i2)))
            currScore = computeCombinedScore(readsToPrint[i]);
        else
            currScore = 0;
        if (currScore>maxScore)
        {
              maxScore = currScore;
              maxIdx = i;
              clear(evenReads);
        }
        if (currScore == maxScore && maxScore > 0)
        {
              if (length(evenReads) == 0)
                    appendValue(evenReads, maxIdx);
              appendValue(evenReads, i);
        }
    }
    if (length(evenReads)>0)
    {
        int randomPrimary;
        srand (time(NULL));
        randomPrimary = rand() % length(evenReads);
        maxIdx = evenReads[randomPrimary];
    }
    //Print "best" alignment to merged bamFile, remember to turn off secondary flag since it might be present(for equally good reads).
    if (!(hasFlagUnmapped(readsToPrint[maxIdx].i1) && hasFlagUnmapped(readsToPrint[maxIdx].i2)))
    {
        if (readsToPrint[maxIdx].i1.qName != "")
        {
            readsToPrint[maxIdx].i1.flag &= ~BAM_FLAG_SECONDARY;
            writeRecord(mergedHaps,readsToPrint[maxIdx].i1);
        }
        if (readsToPrint[maxIdx].i2.qName != "")
        {
            readsToPrint[maxIdx].i1.flag &= ~BAM_FLAG_SECONDARY;
            writeRecord(mergedHaps,readsToPrint[maxIdx].i2);
        }
        erase(readsToPrint, maxIdx);
    }
    for (unsigned i=0; i<length(readsToPrint); ++i)
    {
        int tagValInt = 0;
        BamTagsDict tagsDict1stRead(readsToPrint[i].i1.tags);
        unsigned tagIdx = 0;
        if (readsToPrint[i].i1.qName != "")
        {
            if (!findTagKey(tagIdx, tagsDict1stRead, "XH"))
                std::cerr << "ERROR: No haplotype tag in first read: " << readsToPrint[i].i1.qName << endl;
            else
            {
                if (!extractTagValue(tagValInt, tagsDict1stRead, tagIdx))
                {
                    std::cerr << "ERROR: There was an error extracting XH from tags in: " << readsToPrint[i].i1.qName << endl;
                    if (tagValInt == 0)
                    {
                        if (!(hasFlagUnmapped(readsToPrint[i].i1) && hasFlagUnmapped(readsToPrint[i].i2)))
                        {
                            if (readsToPrint[i].i1.qName != "")
                                writeRecord(equalSecondary,readsToPrint[i].i1);
                            if (readsToPrint[i].i2.qName != "")
                                writeRecord(equalSecondary,readsToPrint[i].i2);
                        }
                    }
                }
            }
        }
    }
    return maxScore;
}

int main(int argc, char const ** argv)
{
    if (argc < 3)
    {
        cerr << "USAGE: " << argv[0] << " hap1.liftedOver.bam hap2.liftedOver.bam [refAligned.bam]\n";
        return 1;
    }
    printStatus("Starting haplotype-merge.");
    //Open bam-files for both haplotypes
    BamFileIn hap1;
    open(hap1, argv[1]);
    BamFileIn hap2;
    open(hap2, argv[2]);
    //Open reference aligned BAM for comparison if provided
    BamFileIn refAl;
    if (argc == 4)
        open(refAl, argv[3]);
    //Declare file, copy header and write header
    BamFileOut mergedHaps("primary.bam", "wb");
    copyHeader(mergedHaps, hap1);
    writeHeader(mergedHaps);
    BamFileOut hap1specific("hap1.specific.bam", "wb");
    copyHeader(hap1specific, hap1);
    writeHeader(hap1specific);
    BamFileOut hap2specific("hap2.specific.bam", "wb");
    copyHeader(hap2specific, hap2);
    writeHeader(hap2specific);
    BamFileOut equalSecondary("equalSecondary.bam", "wb");
    copyHeader(equalSecondary, hap1);
    writeHeader(equalSecondary);
    map<ReadPairInfo, Pair<BamAlignmentRecord> > readPairInfoToReads;
    ReadPairInfo readPair;
    BamAlignmentRecord record, record2, prevRecord, newReadNameRecord, refAlRecord;
    Pair<bool, String<Pair<BamAlignmentRecord> > > allUnalAndReadsToPrint;
    String<Pair<BamAlignmentRecord> > readsForOutput;
    readRecord(prevRecord, hap1);
    unsigned prevCombScore=0, nowCombScore, newBetter=0, oldBetter=0, total=0, nReadsInPrevCombScore = 0;
    addToMap(readPairInfoToReads, prevRecord, 1);
    while (readRecord(record, hap1) && readRecord(record2, hap2))
    {
        while(record.qName == prevRecord.qName)
        {
            addToMap(readPairInfoToReads, record,1);
            if (!readRecord(record, hap1))
                break;
        }
        //Need to process this at end of this loop and set as prevRecord
        newReadNameRecord = record;
        while(record2.qName == prevRecord.qName)
        {
            addToMap(readPairInfoToReads, record2, 2);
            if (!readRecord(record2,hap2))
                break;
        }
        if (readPairInfoToReads.size()>0)
        {
            allUnalAndReadsToPrint = findUniqueAlignments(readPairInfoToReads, hap1specific, hap2specific);
            readsForOutput = allUnalAndReadsToPrint.i2;
            nowCombScore = setPrimarySecondaryFlags(readsForOutput, mergedHaps, equalSecondary);
            if (argc == 5)
            {
                if (!readRecord(refAlRecord, refAl))
                    break;
                while (refAlRecord.qName == prevRecord.qName)
                {
                    if (!hasFlagSecondary(refAlRecord))
                    {
                        ++nReadsInPrevCombScore;
                        prevCombScore += findAlignmentScore(refAlRecord.tags);
                    }
                    if (!readRecord(refAlRecord, refAl))
                        break;
                }
                if (nReadsInPrevCombScore>2)
                    cout << "The prevCombScore contains too many reads" << endl;
                if (prevCombScore > nowCombScore)
                    ++oldBetter;
                if (nowCombScore > prevCombScore)
                    ++newBetter;
                ++total;
                prevCombScore = 0;
                nReadsInPrevCombScore = 0;
                if (!hasFlagSecondary(refAlRecord))
                    prevCombScore += findAlignmentScore(refAlRecord.tags);
            }
            if (allUnalAndReadsToPrint.i1)
            {
                if (readsForOutput[0].i1.qName != "")
                    writeRecord(mergedHaps,readsForOutput[0].i1);
                if (readsForOutput[0].i2.qName != "")
                    writeRecord(mergedHaps,readsForOutput[0].i2);
            }
        }
        readPairInfoToReads.clear();
        addToMap(readPairInfoToReads, newReadNameRecord, 1);
        addToMap(readPairInfoToReads, record2, 2);
        prevRecord = newReadNameRecord;
    }
    std::ostringstream msg;
    msg << "Statistics: ";
    printStatus(msg);
    msg.str("");
    msg.clear();
    msg << "Direct alignment to reference has a higher alignment score in: " << oldBetter << " cases.";
    printStatus(msg);
    msg.str("");
    msg.clear();
    msg << "Alignment to personal genome followed by a liftover to reference genome has a higher alignment score in: " << newBetter << " cases.";
    printStatus(msg);
    msg.str("");
    msg.clear();
    msg << "The total number of cases is: " << total;
    printStatus(msg);
    printStatus("Finished haplotype-merge.");
    return 0;
}
