__author__ = 'armartin (modified by dangliu)'
#takes in rephased alleles rfmix Viterbi file and outputs ASPCA inputs files

import argparse
import os
import gzip
from datetime import datetime
import time

def current_time():
    return '[' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'

def open_file(filename):
    """
    Open files regardless of their gzip status
    """
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file

def check_anc(anc): #TO DO: take in a number corresponding to the ancestral class of interest
    sum_anc = 0
    for i in range(2,len(anc)):
        sum_anc += int(anc[i])
    #print sum_anc
    #print len(anc) - 2
    if sum_anc == 0 or sum_anc == (len(anc) - 2): #check len(keep)
        print str(sum_anc) + ': monomorphic'
        #print anc
        return False
    else:
        return True
    
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def main(args):
    """
    find admixed versus reference panel individuals
    write AS-phased output
    """
    ind_order = [] #reference samples in order to keep
    all_inds = []
    inds = open_file(args.inds)
    classes = open_file(args.classes)
    classes = classes.readline().strip().split()
    classes = classes[0::2] #get odds
    ind_class = {}
    
    keep = set()
    if args.keep is not None:
        filtered = open_file(args.keep)
        for line in filtered:
            keep.add(line.strip())
    
    ## get all admixed inds
    i = 0
    for line in inds:
        line = line.strip()
        all_inds.append(line)
        ind_class[line] = classes[i]
        i += 1

    ind_order = []
    if len(keep) == 0:
        #ind_order = all_inds #fix this to allow for keep file # two ref to same list...
        ind_order = all_inds[:] # copy the original list so that the in_order below won't add up to it
    else:
        for ind in all_inds:
            if ind in keep:
                ind_order.append(ind)
    #print ind_order#
    #print ind_class
    
    ## write beagle header files for admixed and ancestral beagle files
    out_anc = open(args.out + '_anc.beagle', 'w')
    out_adm = open(args.out + '_adm.beagle', 'w')
    out_anc.write('I\tid\t')
    out_adm.write('I\tid\t')
    vit_all = {}
    vit_prob = {}

# Old
#    for ind in all_inds:
#        if ind_class[ind] == '0':
#            out_adm.write(ind + '_A\t' + ind + '_B\t')
#            if len(keep) != 0: # so, the admixed ind number won't be exceed the actual one when there is no "keep" argument
#                ind_order.append(ind) #add admixed samples to ind_order # if not changed above..it will add in all_inds indefinitely unendlich -.-
#        else:
#            if len(keep) == 0:
#                out_anc.write(ind + '_A\t' + ind + '_B\t')
#            elif ind in keep:
#                out_anc.write(ind + '_A\t' + ind + '_B\t')

    for ind in all_inds:
        if ind_class[ind] == '0':
            if len(keep) == 0: # so, the admixed ind number won't be exceed the actual one when there is no "keep" argument
                out_adm.write(ind + '_A\t' + ind + '_B\t')
            elif ind in keep: # make keep can be applied to admixed ind too
                out_adm.write(ind + '_A\t' + ind + '_B\t')
        else:
            if len(keep) == 0:
                out_anc.write(ind + '_A\t' + ind + '_B\t')
            elif ind in keep:
                out_anc.write(ind + '_A\t' + ind + '_B\t')
        vit_all[ind + '_A'] = []
        vit_all[ind + '_B'] = []
    out_adm.write('\n')
    out_anc.write('\n')
    
    #write all markers that are not monomorphic in reference panel
    markers = open_file(args.markers)
    alleles = open_file(args.alleles)
    vit = open_file(args.vit)
    fbk = open_file(args.vit.replace('Viterbi.txt', 'ForwardBackward.txt'))
    out_vit = open(args.out + '.vit', 'w')
    out_markers = open(args.out + '.markers', 'w')
    marker_count = 0
    for line in markers:
        line = line.strip().split()
        markers = alleles.readline().strip().split()
        vit_line = vit.readline().strip().split()
        fbk_line = fbk.readline().strip().split()
        #fbk_chunked = chunker(fbk_line, 3) # for 3 ref -.-
        fbk_chunked = chunker(fbk_line, 2) # for 2 ref
        max_chunks = []
        for chunk in fbk_chunked:
            max_chunks.append(max(map(float, chunk)))
        #print line
        #print markers
        anc = ['M', line[2]]
        adm = ['M', line[2]]
        i=0
        #print markers[0]
        #print vit_line

# Old
#        for ind in all_inds: #this would be a problem because indexing is off
#            if ind_class[ind] == '0':
#                #if ind_class[ind] == '0' and ind in keep:
#                adm.append(markers[0][i])
#                adm.append(markers[0][i+1])
#            elif ind in ind_order:
#                #elif ind.startswith('SA') and ind in keep:
#                anc.append(markers[0][i])
#                anc.append(markers[0][i+1])
#            i += 2

        for ind in all_inds: #this would be a problem because indexing is off
            if len(keep) == 0:
                if ind_class[ind] == '0':
                    adm.append(markers[0][i])
                    adm.append(markers[0][i+1])
                elif ind in ind_order:
                    anc.append(markers[0][i])
                    anc.append(markers[0][i+1])
            else:
                if ind_class[ind] == '0' and ind in ind_order:
                    adm.append(markers[0][i])
                    adm.append(markers[0][i+1])
                elif ind_class[ind] != '0' and ind in ind_order:
                    anc.append(markers[0][i])
                    anc.append(markers[0][i+1])
            i += 2
        i=0
        multimorphic = check_anc(anc)
        if multimorphic: #by default, check all for monomorphic, alternative, check subset for monomorphic
            for ind in all_inds:
                if ind in ind_order:
                    #print ind + ' ' + vit_line[i]
                    if max_chunks[i] > args.fbk_threshold:
                        vit_all[ind + '_A'].append(vit_line[i])
                    else:
                        vit_all[ind + '_A'].append('-9')
                    if max_chunks[i+1] > args.fbk_threshold:
                        vit_all[ind + '_B'].append(vit_line[i+1])
                    else:
                        vit_all[ind + '_B'].append('-9')
                i += 2
            out_adm.write('\t'.join(adm) + '\n')
            out_anc.write('\t'.join(anc) + '\n')
            out_markers.write('window' + str(marker_count) + '\t' + line[2] + '\n')
            marker_count += 1
        
    for ind in ind_order:
        out_vit.write(ind + '_A ' + ' '.join(vit_all[ind + '_A']) + '\n')
        out_vit.write(ind + '_B ' + ' '.join(vit_all[ind + '_B']) + '\n')
    
    out_adm.close()
    out_anc.close()
    out_vit.close()

if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    parser.add_argument('--alleles', help='alleles file produced by RFMix', required=True)
    parser.add_argument('--vit', help='viterbi file produced by RFMix', required=True)
    parser.add_argument('--markers', required=True)
    parser.add_argument('--inds', required=True)
    parser.add_argument('--classes', required=True)
    parser.add_argument('--out', required=True)
    
    parser.add_argument('--fbk_threshold', default=0.99, type=float)
    parser.add_argument('--keep') # Use a list (IID per line) to keep the individuals want to use for aspca
    
    args = parser.parse_args()
    main(args)

# last_v20201203