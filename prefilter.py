"""
Provides various data pre-filters before running PCMA..
"""

import getopt
import sys

def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

                          
    long_options_list = ['in=','out=','filter_weights=','only_zs','only_weights','filter_nan',]

    p_dict = {'in':None,'out':None,'only_zs':False, 'only_weights':False, 'filter_weights':None, 
              'filter_nan':False}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
    
        except:
            print "Some problems with parameters.  Please read the usage documentation carefully."
            print "Use the -h option for usage information."
#             traceback.print_exc()
#             print __doc__
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt == "--in": p_dict['in'] = arg
            elif opt == "--out": p_dict['out'] = arg
            elif opt == "--filter_weights": p_dict['filter_weights'] = float(arg.split(','))
            elif opt == "--filter_nan": p_dict['filter_nan'] = True
            elif opt == "--only_zs": p_dict['only_zs'] = True
            elif opt == "--only_weights": p_dict['only_weights'] = True
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


def only_zs(in_file,out_file):
    with open(in_file,'r') as inf:
        header = inf.next()
        header = header.split()
        weight_i = 0
        while header[weight_i][-7:]!='weights':
            weight_i += 1
        
        with open(out_file,'w') as outf:
            out_str = '\t'.join(header[:weight_i])
            outf.write(out_str+'\n')
            for line in inf:
                l = line.split()
                out_str = '\t'.join(l[:weight_i])
                outf.write(out_str+'\n')
    

if __name__=='__main__':
    p_dict = parse_parameters()
    assert p_dict['in'] is not None, 'Input file is missing.'
    in_file = p_dict['in']
    assert p_dict['out'] is not None, 'Output file is missing.'
    out_file = p_dict['out']

    only_zs(in_file,out_file)
    
    