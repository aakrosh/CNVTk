#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "sam.h"

#include "errors.h"
#include "time.h"
#include "memalloc.h"
#include "files.h"
#include "sequences.h"
#include "hashtable.h"

static bool debug_flag = FALSE;

time_t t0;

typedef int32_t i32;
typedef uint64_t u64;

static void usage()
{
    fprintf(stderr,
    "\nusage:\n"
    "\tcount_in_bins [options] ref.fa ref.mappability alignments.bam\n\n"

    "where the options are:\n"
    "\t-h: print usage and quit\n"
    "\t-d: print debug information\n"
    "\t-b: the number of bases in a single bin [auto]\n\n"

    "Create a file where we report the following:\n"
    "a) chromosome\n"
    "b) interval start\n"
    "c) interval end\n"
    "d) number of sequences that start in the bin\n"
    "e) number of bases in bin that are in [ACGT]\n"
    "f) GC content of the bin\n\n"

    "for non-overlapping windows of size specified by -b mappable bases. If\n"
    "-b is not specified we select the size automatically by requiring that\n" 
    "at least 100 fragments on average cover a bin composed of mappable bases\n" 
    "The mappable segments are constructed from a mappability file, that \n"
    "should be produced using gem-mappability\n"
    "(http://algorithms.cnag.cat/wiki/The_GEM_library#Documentation). The\n"
    "manuscript that discusses the mappability is here:\n"
    "\"http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377.\"\n\n)");
}

typedef struct chrcoverage_t
{
    u64 length;
    uchar* map; // taking one byte per base, this can be reduced further. $$$
    uchar* cov; // this allows me to store a coverage upto 255 per base
    uchar* seq;
} chrcoverage;


static hashtable* ReadReference(const char* const refname)
{
    hashtable* reference = new_hashtable(12);

    sequence* sp = read_fasta_sequence(refname);
    
    while(sp != NULL){
        // allocate a coverage array for the sequence
        chrcoverage* cov = ckallocz(sizeof(chrcoverage));
        cov->length   = strlen((char*)sp->sequence);
        cov->map = ckallocz(strlen((char*)sp->sequence));
        cov->cov = ckallocz(strlen((char*)sp->sequence));
        cov->seq = ckallocz(strlen((char*)sp->sequence)+1);
        memcpy(cov->seq, sp->sequence, cov->length);

        // if the name of the sequence has more than one tokens, just use the
        // first token in the name
        int i = 0;
        while((sp->header[i] != '\n') && 
              (sp->header[i] != 0)    && 
              (sp->header[i] != '\t') && 
              (sp->header[i] != 32)) i++;
        sp->header[i] = 0;

        add_hashtable(reference,(char*)sp->header,strlen((char*)sp->header),cov);
        sp = get_next_sequence(sp);
    } 

    return reference;
}

static void CreateRefMap(const char* const refName,
                         const char* const mapName,
                         hashtable* const reference)
{
    /* Return an object that lists the unique regions of the reference.
 
        The format of the input file includes description of the encoding used,
        followed by fasta like format.   
        ~~ENCODING
        ' '~[0-0]
        '!'~[1-1]
        '"'~[2-2]
         ...
        '&'~[6-7]
         ...
        ~chr17
        !!!!!!!!!!!!!!!!!!!!!!!!!""!!!!!!!!!!%%*&&00,00/66/140.-,,
        ,/04:41237?CDDDD??>;7*$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!%$$#!!"#$&%"#.7>:><9898899999:.3---1111"!!!!##!!
        ...
    */
    
    FILE* fp = ckopen(mapName, "r");
    size_t n = 1;
    char* line = ckallocz(n*sizeof(char));
    bool encoding_flag = FALSE;
    bool map_flag = FALSE;
    char sentinel = '~';

    chrcoverage* cc;
    int index = 0;
    int num_mappable = 0, num_total = 0;
    signed long ret;    

    while ((ret = getline(&line, &n, fp)) != -1) {
        if  (line[0] == '~') {
            // beginning of a subsection
            if (strncmp(line, "~~ENCODING", 10) == 0) {
                encoding_flag = TRUE;  
            } else if (map_flag == TRUE) {
                cc = (chrcoverage*)must_find_hashtable(reference, line+1, strlen(line+1)-1);
                index = 0;
            }
        } else {
            if (encoding_flag == TRUE) {
                if (strncmp(line+3, "~[1-1]",6) == 0) {
                    sentinel = line[1];
                    encoding_flag = FALSE;
                    map_flag = TRUE;
                }
            } else {
                int i, j;
                for (j =0; j < ret; j++) {
                    if (line[j] == '\n') continue;
                    if (line[j] == sentinel) {
                        cc->map[index] = '1';  
                        num_mappable += 1;
                    }
                    index += 1;
                    num_total += 1;
                }
            }
        }
    }

    ckfree(line);
    fclose(fp);

    fprintf(stderr, "%d (%2.2f%%) mappable bases in the genome.\n", 
        num_mappable, num_mappable * 100.0 / num_total);

    //bin* iter;
    //bin* next;
    //
    //u64 sum = 0;
    //u64 num = 0;
    //for(int i = 0; i < reference->size; i++){
    //    iter = reference->bins[i];
    //    while(iter){
    //        next = iter->next;

    //        chrcoverage* chrcov = (chrcoverage*)iter->val;
    //        for (u64 j = 0; j < chrcov->length; j++) {
    //            if (chrcov->map[j] == '1') {
    //                printf("%s\t%"PRIu64"\n", iter->name, j);
    //            }
    //        }

    //        iter = next;
    //    }
    //}   
   
}

static int CreateCoverageMap(const char* const refName,
                             const char* const bamName,
                             hashtable* const reference)
{
    int status;
    samFile *in = sam_open(bamName, "r");
    bam_hdr_t *hdr = NULL;
    if (in == NULL) {
        status = EXIT_FAILURE;
        return 0;
    }
    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        status = EXIT_FAILURE;
        goto clean;
    }

    int ret;
    bam1_t *b = bam_init1();
    u64 numread = 0; // number of reads analyzed

    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        numread += 1;
        if ((numread % 10000000) == 0) {
            fprintf(stderr, "Processed %"PRIu64" reads\n", numread);
        }
        if (1 == debug_flag) {
            fprintf(stderr, "Read name : %s\n", bam_get_qname(b));
        }

        // ignore if this is a zero length read (have seen it in some cases
        // where the reads were clipped by another tool. also ignore all
        // secondary or supplementary or QC failed alignments for now. 
        // ignore unmapped reads as well
        if (b->core.l_qseq == 0) continue;
        if (((b->core.flag & 0x4) == 0x4) || 
            ((b->core.flag & 0x100) == 0x100) || 
            ((b->core.flag & 0x200) == 0x200) || 
            ((b->core.flag & 0x400) == 0x400) ||
            ((b->core.flag & 0x800) == 0x800)){
            continue;
        }

        // if this is paired, then I register one vote for the fragment.
        if ((b->core.flag & 0x1) == 0x1) {
            if ((b->core.flag & 0x40) == 0x40) {
                chrcoverage* cov = must_find_hashtable(reference,
                               hdr->target_name[b->core.tid],
                               strlen(hdr->target_name[b->core.tid]));
                cov->cov[b->core.pos] += 1; 
                if(cov->cov[b->core.pos] == 251) cov->cov[b->core.pos] = 250;
            }
        } else {
            chrcoverage* cov = must_find_hashtable(reference,
                               hdr->target_name[b->core.tid],
                               strlen(hdr->target_name[b->core.tid]));
            cov->cov[b->core.pos] += 1;
            if(cov->cov[b->core.pos] == 251) cov->cov[b->core.pos] = 250;
        }
    }

clean:
    if (hdr != NULL) bam_hdr_destroy(hdr);
    if (hts_close(in) != 0)
        status = EXIT_FAILURE;

    //bin* iter;
    //bin* next;
    //
    //u64 sum = 0;
    //u64 num = 0;
    //for(int i = 0; i < reference->size; i++){
    //    iter = reference->bins[i];
    //    while(iter){
    //        next = iter->next;

    //        chrcoverage* chrcov = (chrcoverage*)iter->val;
    //        for (u64 j = 0; j < chrcov->length; j++) {
    //            if (chrcov->cov[j] > 0) {
    //                printf("%s\t%"PRIu64"\t%d\n", iter->name, j, chrcov->cov[j]);
    //            }
    //        }

    //        iter = next;
    //    }
    //}   
 

    return status;
}

static float CalculateAverageCoverage(hashtable* const reference, const int bs)
{
    float average_coverage = 0;

    bin* iter;
    bin* next;
    
    u64 sum = 0;
    u64 num = 0;
    for(int i = 0; i < reference->size; i++){
        iter = reference->bins[i];
        while(iter){
            next = iter->next;

            int binSum = 0; 
            int binIndex = 0;

            chrcoverage* chrcov = (chrcoverage*)iter->val;
            for (u64 j = 0; j < chrcov->length; j++) {
                if (chrcov->map[j] == '1') {
                    binSum += chrcov->cov[j];
                    binIndex += 1;
                    if (binIndex == 100) {
                        // only consider non-zero bins
                        if (binSum > 0) {
                            sum += binSum;
                            num += 1;
                        }
                        binIndex = 0;
                        binSum = 0;
                    }
                }
            }

            iter = next;
        }
    }   

    
    average_coverage = sum * 1.0 / num;
    return average_coverage;
}

static int CalculateAGoodBinSize(const char* const refName,
                                 hashtable* const reference)
{
    int binSize = 100;

    float average_cov = CalculateAverageCoverage(reference, binSize);
    fprintf(stderr, "Average coverage with 100 mappable bases: %2.2f\n", average_cov);    

    if (average_cov >= 100)
        binSize = 100;
    else
        binSize = 100 * (100 * 1.0 / average_cov);

    return binSize;
}

static int CalculateGC(const chrcoverage* const chrcov,
                       const u64 start,
                       const u64 end)
{
    int a = 0, c = 0, g = 0, t = 0;
    for (u64 i = start; i < end; i++) {
        if (chrcov->map[i] == '1') {
            switch(toupper(chrcov->seq[i])) {
                case 'A':
                    a += 1;
                    break;
                case 'C':
                    c += 1;
                    break;
                case 'G':
                    g += 1;
                    break;
                case 'T':
                    t += 1;
                    break;
                default:
                    break;
            }
        } 
    }

    return (g+c) * 100.0 / (a+c+g+t);
}

static void PrintBinInfo(const char* const refName,
                         hashtable* const reference,
                         const int binSize)
{
    bin* iter;
    bin* next;
    int gc;
    
    for(int i = 0; i < reference->size; i++){
        iter = reference->bins[i];
        while(iter){
            next = iter->next;

            int binSum = 0; 
            int binIndex = 0;
            u64 start = 0;

            chrcoverage* chrcov = (chrcoverage*)iter->val;

            for (u64 j = 0; j < chrcov->length; j++) {
                if (chrcov->map[j] == '1') {

                    binSum += chrcov->cov[j];
                    binIndex += 1;
                    if (binIndex == binSize) {
                        gc = CalculateGC(chrcov, start, j+1);
                        printf("%s\t%"PRIu64"\t%"PRIu64"\t%d\t%d\t%d\n", 
                            iter->name, start, j+1, binSum, binSize, gc);
                        binIndex = 0;
                        binSum = 0;
                        start = j + 1;
                    }
                }
            }
    
            if (binIndex > 0) {
                gc = CalculateGC(chrcov, start, chrcov->length);
                printf("%s\t%"PRIu64"\t%"PRIu64"\t%d\t%d\t%d\n", 
                    iter->name, start, chrcov->length, binSum, binIndex, gc);
            }

            iter = next;
        }
    }   
}

static void count_in_bins(const char* const refName,
                          const char* const mapName,
                          const char* const bamName,
                          int binSize)
{
    // lets read the reference sequence. The reference sequence will be a
    // hashtable where the key is the name of the chromosome and the value will
    // be a char array of size equal to  the length of the chromosome sequence.
    hashtable* reference = ReadReference(refName);
    timestamp("Read the reference sequence.\n");

    // create a map of bases that are  uniquely mappable
    CreateRefMap(refName, mapName, reference);
    timestamp("Created a map of unique locations in the genome.\n");

    // create a coverage map of all locations from the BAM file
    CreateCoverageMap(refName, bamName, reference);
    timestamp("Read the coverage information from the BAM file.\n");

    // calculate a bin size that will ensure that at least 100 fragments cover a
    // bin on average
    if (binSize == 0) 
        binSize = CalculateAGoodBinSize(refName, reference);
    fprintf(stderr, "Using a bin size of %d mappable bases.\n", binSize);

    // using that bin size, print out the information about the bins
    PrintBinInfo(refName, reference, binSize);
    timestamp("Done with the counting.\n");
}

int main (int argc, char **argv)
{
    int binsize = 0;
    argv0 = "count_in_bins";    

    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "hdb:")) != -1)
        switch (c) {
            case 'h':   
                usage();
                return EXIT_SUCCESS;
            case 'd':
                debug_flag = TRUE;
                break;
            case 'b':
                binsize = atoi(optarg);
                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", 
                    optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", 
                    optopt);
                return EXIT_FAILURE;
            default:
                abort ();
        }
    
    if ((argc - optind) != 3) {
        usage();
        return EXIT_FAILURE;
    }

    t0 = time(0);
    count_in_bins(argv[optind], argv[optind+1], argv[optind+2], binsize);
    print_usage();
    
    return EXIT_SUCCESS;
}
