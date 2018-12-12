#include <stdlib.h>
#include <map>
#include <algorithm>
#include <seqan/bam_io.h>
#include <gapped_sequence.h>

using namespace std;
using namespace seqan;

void printStatus(const char * message)
{
        // Get the current date and time.
        char timestamp[80];
        time_t now = time(0);
        struct tm tstruct;
        tstruct = *localtime(&now);
        strftime(timestamp, sizeof(timestamp), "[bamLiftOver %Y-%m-%d %X] ", &tstruct);

        // Print time and message.
        std::cerr << timestamp << message << std::endl;
}

void printStatus(std::ostringstream & message)
{
        std::string msg = message.str();
        printStatus(toCString(msg));
}

int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        cerr << "USAGE: " << argv[0] << " readsToPnRef.bam pnRefToRealRef.chain outputFile.bam\n";
        return 1;
    }

    //Read chainfile into a map.
    CharString chainfile = argv[2];
    map<CharString,Pair<GappedSequence> > chrToChain;
    readChainFile(chrToChain, chainfile);
    printStatus("Finished reading chain-file.");

    //Open bamfiles
    HtsFile readsToPnRef(argv[1], "r");
    HtsFile bamFileOut(argv[3], "wb");

    //Copy and update header
    copyHeader(bamFileOut, readsToPnRef);
    vector<CharString> chromNames;
    updateHeader(bamFileOut.hdr, chrToChain, chromNames);

    //Write updated header to outputFile.
    writeHeader(bamFileOut);
    printStatus("Finished updating and printing header.");

    BamAlignmentRecord record;
    unsigned oldStartPos;

    while (readRecord(record, readsToPnRef))
    {
        //Read bamfile: If read is unaligned or its chr is not in the chain file, do nothing and write it to output.
        //Might have to update mate position even though read is unaligned
        if (hasFlagMultiple(record))
            updateMatePosition(record, chrToChain, chromNames);
        if (record.rID == -1)
            writeRecord(bamFileOut,record);
        else
        {
            oldStartPos = record.beginPos;
            //Update coordinates and if record.beginPos is updated, then update cigarString
            if (updateBamPosition(record, chrToChain, chromNames))
            {
                 updateCigarString(record, chrToChain, chromNames, oldStartPos);
            }
            //Print updated record to output
            writeRecord(bamFileOut, record);
        }
    }
    printStatus("Finished lift-over.");
    return 0;
}
