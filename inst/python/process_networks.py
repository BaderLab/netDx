#!/usr/bin/env python

import sys
import os

from os.path import join

batch_file_name = sys.argv[1]

organisms = {}
organism_ids = {}
organism_data = {}
organisms_by_group = {}

# GeneMania db files created by this script.
genes_file = open('GENES.txt', 'w')
nodes_file = open('NODES.txt', 'w')
gene_data_file = open('GENE_DATA.txt', 'w')
group_file = open('NETWORK_GROUPS.txt', 'w')
colour_file = open('colours.txt', 'w')
network_file = open('NETWORKS.txt', 'w')
metadata_file = open('NETWORK_METADATA.txt', 'w')
organism_file = open('ORGANISMS.txt', 'w')

# dummy variable initialization
organism_id = 1
naming_source_id = 1

group_id = 1
groups = {}

network_id = 1

try:
    os.mkdir('INTERACTIONS')
except:
    pass

def handle_ids(naming_source_id, organism_id, id_file_name):
    ids = []
    id_map = {}
    for line in open(id_file_name, 'rU'):
        ids.append(line.strip())
    for index, symbol in zip(range(len(ids)), ids):
        node_id = index + 1
        print >> nodes_file, '%(node_id)d\t%(symbol)s\t%(node_id)d\t%(organism_id)d' % locals()
        print >> genes_file, '%(node_id)d\t%(symbol)s\tN/A\t%(naming_source_id)d\t%(node_id)d\t%(organism_id)d\t0' % locals()
        print >> gene_data_file, '%(node_id)d\t%(symbol)s\t\t' % locals()
        id_map[symbol] = node_id
    return id_map

def handle_synonyms(organism_id, id_map):
    ids = set(id_map.values())
    synonym_file_name = '%d.synonyms' % organism_id
    synonym_file = open(synonym_file_name, 'w')
    try:
        for id in ids:

            print >> synonym_file, '%d\t%d' % (id, id)
    finally:
        synonym_file.close()
    
def handle_organism(data):
    global organism_id
    global naming_source_id
    ids = handle_ids(naming_source_id, organism_id, data[1])
    handle_synonyms(organism_id, ids)
    
    organism_ids[organism_id] = ids
    organism_data[organism_id] = data
    organisms[data[2]] = organism_id
    organism = [str(organism_id), data[2], data[3], data[4], '-1', data[5]]
    print >> organism_file, '\t'.join(organism)
    organism_id += 1
    
def handle_group(data):
    global group_id
    global group_file
    global colour_file
    global groups
    group = [str(group_id), data[1], data[2], data[3], str(organisms[data[5]])]
    print >> group_file, '\t'.join(group)
    print >> colour_file, '%s\t%s' % (data[2], data[4])
    groups[data[2]] = group_id
    organisms_by_group[group_id] = int(organisms[data[5]])
    group_id += 1

def count_interactions(organism_id, file_name):
    symbols = organism_ids[organism_id]
    interactions = 0
    
    if file_name.endswith('.profile'):
        return 0
    
    for line in open(file_name, 'rU'):
        data = line.strip().split('\t')
        if data[0] in symbols and data[1] in symbols:
            interactions += 1
    return interactions

def write_network(organism_id, file_name, network_id):
    symbols = organism_ids[organism_id]
    id_indexes = [0, 1]# SP modified this for situations
						   # where pairwise interaction networks are 
						   # provided
    if file_name.endswith('.profile'):
        id_indexes = [0]
        output_file_name = '%d.%d.profile' % (organism_id, network_id)
        output_path = 'profiles'
    else:
        output_file_name = '%d.%d.txt' % (organism_id, network_id)
        output_path = 'INTERACTIONS'
        
    output_file = open(join(output_path, output_file_name), 'w')
    for line in open(file_name, 'rU'):
        data = line.strip().split('\t')
#	print data[0];
#	print file_name;
	#print symbols;

        for i in id_indexes:             
            data[i] = str(symbols[data[i]])
        print >> output_file, '\t'.join(data)        
    
def handle_network(data):
    global network_id
    global network_file
    global metadata_file
    global groups
    group_id = groups[data[4]]
    organism_id = organisms_by_group[group_id]
    network = [str(network_id), data[2], str(network_id), data[3], '0', str(group_id)]
    print >> network_file, '\t'.join(network)
    
    file_name = data[1]
    interactions = count_interactions(organism_id, file_name)
    
    metadata = [str(network_id), '', '', '', '', '', '', '', '', '', str(interactions), '', '', '0', '', '', '', '', '']
    print >> metadata_file, '\t'.join(metadata)
    write_network(organism_id, file_name, network_id)
    network_id += 1

def write_configuration(file_name, organism_data):
    output_file = open(file_name, 'w')
    print >> output_file, '''[FileLocations]
generic_db_dir = .'''

    organisms = ['org_%d' % id for id in organism_data.keys()]
    print >> output_file, '''[Organisms]
organisms = %s''' % (','.join(organisms),)

    for id, data in organism_data.items():
        print >> output_file, '''[org_%d]
gm_organism_id = %d
short_name = %s
common_name = %s
''' % (id, id, data[2], data[3])

def write_naming_sources():
    output_file = open('GENE_NAMING_SOURCES.txt', 'w')
    print >> output_file, '''1	Other	1	Other
2	Entrez Gene ID	0	Entrez Gene ID'''

for line in open(batch_file_name, 'rU'):
    data = line.strip().split('\t')
    if data[0] == 'organism':
        handle_organism(data)
    if data[0] == 'group':
        handle_group(data)
    if data[0] == 'network':
        handle_network(data)

write_naming_sources()
write_configuration('db.cfg', organism_data)

