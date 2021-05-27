#!/usr/bin/env python3
#import pandas as pd
import sys, getopt, re, yaml

def get_ichorPath(rlib_path, config):
    #print(str(rlib_path))
    file=open(str(rlib_path), mode='r',newline="\n")
    rlib_path = file.read()
    #print(str(rlib_path))
    extdata  = str(rlib_path).rstrip() + "/ichorCNA/extdata/"
    
    # Setup centromere file name (e.g. GRCh37 instead of hg19)
    if config['common']['build'] == 'hg19':
        cen_file = 'GRCh37.p13_centromere_UCSC-gapTable.txt'
    elif config['common']['build'] == 'hg38':
        cen_file = 'GRCh38.GCA_000001405.2_centromere_acen.txt'
        #cen_file = 'cytoBand_hg38'
    
    # Get the 500kb or 1Mb window annotation
    window_size = config['params']['readcounter']['window']
    window_size = int(window_size / 1000)       # size in kb
    if window_size >= 1000:
        window_size_simple = str(int(window_size / 1000)) + "Mb"
    else:
        window_size_simple = str(window_size) + "kb"
    
    # assemble wig file prefix
    # Expected Format: [gc/map]_[hg19/hg38]_[window_size]kb.wig
    wig_file = "_" + config['common']['build'] + "_" + str(window_size) + "kb.wig"
    normal_file = "HD_ULP_PoN_" + window_size_simple + "_median_normAutosome_mapScoreFiltered_median.rds"
    map_path    = extdata + "map" + wig_file
    gc_path     = extdata + "gc" + wig_file
    cen_path    = extdata + cen_file
    normal_path = extdata + normal_file
    return { "map":map_path, "gc":gc_path, "cen":cen_path, "norm":normal_path }

def get_ichorChrs(chr_path, config):
    #print(str(chr_path))
    file=open(str(chr_path), mode='r',newline="\n")
    chrs = file.read()
    
    chrs        = re.sub("chr", "", chrs.rstrip())
    chr_train   = "c(" + re.sub(",X.*$", "", chrs) + ")"
#    chrs        = "c(" + re.sub(",Y.*$", "", chrs) + ")"
    chrs        = "c(" + re.sub(",X.*$", "", chrs) + ")"
    
    return {"all":chrs, "train":chr_train}

def main(argv):
    # Config file
    configfile= "config/config.yaml"
    file=open(configfile)
    config=yaml.load(file, Loader=yaml.FullLoader)
    
    inputfile = ''
    returnarg = ''
    returnval = ''
    try:
        opts, args = getopt.getopt(argv, "hi:r:",["ifile=", "returnarg="])
    except getopt.GetoptError:
        print('parse_paths.py -i <inputfile> -r <returnarg>')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print('parse_paths.py -i <inputfile> -r <returnarg>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-r", "--returnarg"):
            returnarg = arg
    
    if returnarg in ['map', 'gc', 'cen', 'norm']:
        returnval = get_ichorPath(inputfile, config)[returnarg]
    elif returnarg in ['all', 'train']:
        returnval = get_ichorChrs(inputfile, config)[returnarg]
    print(returnval)

if __name__ == "__main__":
    main(sys.argv[1:])
