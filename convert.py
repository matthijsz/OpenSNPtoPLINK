import _csv
import os
import pandas as pd
import snpy as sn
from plinkfiles import PlinkFiles

compress = False  # Set this to True if you run into memory errors, though that will significantly slow things down
# Setting this to True will compress the binary genotype data while not in use
# The main issue is for each SNP (per participant) it will then have to decompress and recompress the genotypes,
#  significantly slowing everything down.

autosomes = list(range(1, 23))
fam, logdata, file_queue = {}, {}, []
n_participants = 0
all_files = os.listdir('RAW')
# First assess (a) which files will run succesfully, and (b) how many participants there are
# In the meanwhile we can also generate the .fam file
for nfile, filename in enumerate(all_files):
    if filename.startswith('user'):
        try:
            # Trying to catch file errors
            snps = sn.parse("RAW/{}".format(filename))
            for row in snps:
                pass
        except UnicodeError:
            # Some files produce UnicodeDecodeError, probably different format? Skip those
            print('{}: UnicodeError'.format(filename))
            logdata[filename] = 'UnicodeError'
            continue
        except TypeError:
            # No idea what happens here, but skip em
            print('{}: TypeError'.format(filename))
            logdata[filename] = 'TypeError'
            continue
        except ValueError:
            # No idea what happens here, but skip em
            print('{}: ValueError'.format(filename))
            logdata[filename] = 'ValueError'
            continue
        except RuntimeError:
            # PyVCF failed to launch
            print('{}: PyVCF failed'.format(filename))
            logdata[filename] = 'PyVCFError'
            continue
        except IndexError:
            # Empty files ???
            print('{}: IndexError'.format(filename))
            logdata[filename] = 'IndexError'
            continue
        except _csv.Error:
            # How much more of this?
            print('{}: CSVError'.format(filename))
            logdata[filename] = 'CSVError'
            continue
        famdata = filename.split('_')
        if famdata[0].replace('user', '') not in fam.keys():
            n_participants += 1
            print('Verifying file {}/{} - subject {}: {}'.format(nfile, len(all_files), n_participants, filename))
            if famdata[-1].split('.') == 'unkown':
                sex = 0
            else:
                sex = 1 if famdata[-1].split('.') == 'XX' else 2
            fam[famdata[0].replace('user', '')] = {
                'FID': famdata[0].replace('user', ''),
                'IID': famdata[0].replace('user', ''),
                'PID': 0,
                'MID': 0,
                'sex': sex,
                'phen': 0
            }
            file_queue.append(filename)

# Save .fam file
pd.DataFrame.from_dict(fam, orient='index').to_csv("opensnp.fam", sep='\t', index=False, header=False)

plinkfiles = PlinkFiles(n_participants, compress=compress)
triallelic = []
for n_participant, filename in enumerate(file_queue):
    print('Converting file {}/{}: {}'.format(n_participant, n_participants, filename))
    n_snps_file = 0
    snps = sn.parse("RAW/{}".format(filename))
    processed_variantes = {i: {} for i in autosomes}
    for row in snps:
        try:
            chrom = int(row.chromosome)
        except ValueError:
            continue
        if (chrom in autosomes) and (row.name not in triallelic):
            if row.name in processed_variantes[chrom].keys():
                triallelic.append(row.name)
            else:
                plinkfiles.add(n_participant, row.name, chrom, row.position, row.genotype)
                processed_variantes[chrom][row.name] = True
                n_snps_file += 1
        else:
            break  # Assuming all files are stored in order we can stop once we hit a non-autosomal chromosome
    logdata[filename] = 'Succes,{}'.format(n_snps_file)

plinkfiles.remove(triallelic)
# Write .bed and .bim files
plinkfiles.save('opensnp')

# Write log
with open("opensnp.log", 'w') as f:
    f.write("Total subjects: {}\n".format(n_participant))
    f.write("Total SNPs    : {}\n\n".format(len(plinkfiles.bim)))
    succes = {k: v for k, v in logdata.items() if v.startswith('Succes')}
    failed = {k: v for k, v in logdata.items() if not v.startswith('Succes')}
    f.write('FAILED FILES:')
    for k, v in failed.items():
        f.write('\n  {} - {}'.format(k, v))
    f.write('\n\nSUCCSESFULLY PARSED FILES:')
    for k, v in succes.items():
        f.write('\n  {}: extracted {} SNPs'.format(k, v))




