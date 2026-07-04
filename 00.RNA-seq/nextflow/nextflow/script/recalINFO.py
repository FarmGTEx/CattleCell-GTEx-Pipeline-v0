# -*- coding: utf-8 -*-
"""
Created on Wed May 23 17:00:34 2023

@author: Gongmian
"""

import click
import os


@click.command()
@click.option('--invcf', help='input vcf.gz/bcf file')

def main(invcf):
    for line in os.popen(f"bcftools view {invcf}", 'r'):
        if '#' in line:
            print(line.strip())
        else:
            tline = line.strip().split()
            RAF = tline[7].split(';')[0]
            dsind = tline[8].split(':').index('DS')
            gpind = tline[8].split(':').index('GP')

            ds_sum = ds2_sum = ds4_sum = 0

            for samid in tline[9:]:
                samidl = samid.split(':')
                ds = float(samidl[dsind])
                gp0, gp1, gp2 = samidl[gpind].split(',')
                gp1 = float(gp1)
                gp2 = float(gp2)
                
                ds_sum += ds
                ds2_sum += ds * ds
                ds4_sum += gp1 + 4.0 * gp2

            num = len(tline[9:])
            AF = ds_sum / (2 * num)
            if AF == 0 or AF == 1:
                INFO = 1
            else:
                INFO = 1.0 - (ds4_sum - ds2_sum) / (2 * num * AF * (1.0 - AF))
                #INFO = 1.0 - (ds4_sum - ds2_sum) / (2 * (num + refnum) * AF * (1.0 - AF))
            if INFO < 0:
                INFO = 1
            tline[7] = f'{RAF};AF={AF:.6f};INFO={INFO:.6f}'
            
            print('\t'.join(tline))
        

if __name__ == '__main__':
    main()
