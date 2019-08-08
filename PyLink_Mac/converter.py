#!/usr/bin/python
import sys
import numpy as np
import argparse
import re
import itertools
import os
import glob
import time
from math import sqrt

################ CONSTANTS ################
date = "07.05.2018"
global amino_acids, nucleotides, names, WandaLassoProgram, min_loop_length, artifact_gap_size, bridge_max_length, close_residues_cutoff
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'BTC', 'FCY', 'GGL']
nucleotides = ['DG', 'DT', 'DU', 'DA', 'DC','G','T','A','C']
ions = ['CA', 'FE', 'ZN', 'MG']
names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnoprstuvwxyz"
WandaLassoProgram = 'surfacesMyOrient'          # name of the program for lasso analysis
min_loop_length = 5                             # smallest loop size allowed
artifact_gap_size = 6                           # smallest gap size to rise warning and classify structure to artifacts
bridge_max_length = 12                          # maximal length of the bridge to exist
close_residues_cutoff = 4.6                     # which residues treat as close enough

################ FUNCTIONS ################
def dfs(graph, start, end):
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if path and state == end:
            yield '-'.join([start] + path)
            continue
        for next_state in graph[state]:
            if next_state[1] in path: 
                continue
            fringe.append((next_state[1], path+next_state))
def print_namespace(args):              #format namespaces for debug option
    print color.BOLD + color.UNDERLINE + "Passed arguments:" + color.END
    print color.BOLD + "input_file\t" + color.END + args.input_file
    if args.traj : print color.BOLD + "traj\t\t" + color.END + color.GREEN + str(args.traj) + color.END
    if not args.traj : print color.BOLD + "traj\t\t" + color.END + color.RED + str(args.traj) + color.END
    if args.fourcolumn : print color.BOLD + "fourcolumn\t" + color.END + color.GREEN + str(args.fourcolumn) + color.END
    if not args.fourcolumn : print color.BOLD + "fourcolumn\t" + color.END + color.RED + str(args.fourcolumn) + color.END    
    if args.atoms : print color.BOLD + "atoms\t\t" + color.END + str(args.atoms)
    if not args.atoms : print color.BOLD + "atoms\t\t" + color.END + color.GREEN + "CA, C3'" + color.END
    if args.permitted : print color.BOLD + "permitted\t" + color.END + str(args.permitted)
    if not args.permitted : print color.BOLD + "permitted\t" + color.END + color.GREEN + "All" + color.END
    print color.BOLD + "debug\t\t" + color.END + str(args.debug)
    if args.bridges : print color.BOLD + "bridges\t\t" + color.END + color.GREEN + str(args.bridges) + color.END
    if not args.bridges : print color.BOLD + "bridges\t\t" + color.END + color.RED + str(args.bridges) + color.END    
    print color.BOLD + "commands\t" + color.END + str(args.commands)
    if args.extended : print color.BOLD + "extended\t" + color.END + color.GREEN + str(args.extended) + color.END
    if not args.extended : print color.BOLD + "extended\t" + color.END + color.RED + str(args.extended) + color.END    
    if args.sbridge : print color.BOLD + "sbridge\t\t" + color.END + color.GREEN + str(args.sbridge) + color.END
    if not args.sbridge : print color.BOLD + "sbridge\t\t" + color.END + color.RED + str(args.sbridge) + color.END    
    if args.schain : print color.BOLD + "schain\t\t" + color.END + color.GREEN + str(args.schain) + color.END
    if not args.schain : print color.BOLD + "schain\t\t" + color.END + color.RED + str(args.schain) + color.END    
    if args.macro : print color.BOLD + "macro\t\t" + color.END + str(args.macro)
    if not args.macro : print color.BOLD + "macro\t\t" + color.END + color.RED + str(args.macro) + color.END    
    if args.xyz : print color.BOLD + "xyz\t\t" + color.END + str(args.xyz)
    if not args.xyz : print color.BOLD + "xyz\t\t" + color.END + color.RED + str(args.xyz) + color.END  
    print color.BOLD + "alt location\t" + color.END + str(args.altloc)
    print "________________________________\n"

def read_header(pdb_file,args):                          # reading header of PDB file
    input_file = args['input_file']
    debug = args['debug']
    if input_file == 'ssss' : 
        snake()
        sys.exit()
    if not os.path.isfile(input_file):
        print "Cannot find " + input_file + "! Exiting."
        sys.exit()
    for f in glob.glob(input_file + '_*.pdb'): 
        if f[:6] != 'macro_': os.remove(f)            #removing files on start
    for f in glob.glob(input_file + '_*.xyz'): 
        if f[:6] != 'macro_': os.remove(f)            #removing files on start
    with open(input_file, 'r') as f:
        if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) + " File " + input_file + " opened. Analyzing the header."
        for line in f:
            if line[0:6] == 'HEADER': pdb_file.add_header(line)
            if line[0:6] == 'OBSLTE': pdb_file.add_obsolete(line)
            if line[0:6] == 'TITLE ': pdb_file.add_title(line)
            if line[0:6] == 'SPLIT ': pdb_file.add_split(line)
            if line[0:6] == 'CAVEAT': pdb_file.add_caveat(line)
            if line[0:6] == 'COMPND': pdb_file.add_compound(line)
            if line[0:6] == 'SOURCE': pdb_file.add_source(line)
            if line[0:6] == 'KEYWDS': pdb_file.add_keywords(line)
            if line[0:6] == 'EXPDTA': pdb_file.add_expmethod(line)
            if line[0:6] == 'NUMMDL': pdb_file.add_nummdl(line)
        #    if line[0:6] == 'MDLTYP': pdb_file.add_nummdl(line) 
        #    if line[0:6] == 'AUTHOR': pdb_file.add_author(line)
        #    if line[0:6] == 'REVDAT': pdb_file.add_revdat(line)
        #    if line[0:6] == 'SPRSDE': pdb_file.add_sprsde(line)
        #    if line[0:6] == 'JRNL  ': pdb_file.add_journal(line)
            if line[0:10] == 'REMARK 465': pdb_file.add_missing_residue(line)
            if line[0:10] == 'REMARK 470': pdb_file.add_zero_occupancy_residue(line)
        #    if line[0:10] == 'REMARK 475': pdb_file.add_missing_atom(line)
            if line[0:10] == 'REMARK 480': pdb_file.add_zero_occupancy_atom(line)
        #    if line[0:6] == 'DBREF ': pdb_file.add_dbref(line)
        #    if line[0:6] == 'SEQADV': pdb_file.add_seqadv(line)
            if line[0:6] == 'SEQRES': pdb_file.add_residue(line)
        #    if line[0:6] == 'MODRES': pdb_file.add_modres(line)
        #    if line[0:6] == 'HET   ': pdb_file.add_het(line)
        #    if line[0:6] == 'HETNAM': pdb_file.add_hetnam(line)
        #    if line[0:6] == 'HETSYN': pdb_file.add_hetsyn(line)
        #    if line[0:6] == 'FORMUL': pdb_file.add_formul(line)
        #    if line[0:6] == 'HELIX ': pdb_file.add_helix(line)
        #    if line[0:6] == 'SHEET ': pdb_file.add_sheet(line)
            if line[0:6] == 'SSBOND': pdb_file.add_bond(line)
            if line[0:6] == 'LINK  ': pdb_file.add_bond(line)
        #    if line[0:6] == 'CISPEP': pdb_file.add_cispep(line)
        #    if line[0:6] == 'SITE  ': pdb_file.add_site(line)
        #    if line[0:6] == 'CRYST1': pdb_file.add_cryst1(line)
        #    if line[0:6] == 'ORIGX1': pdb_file.add_origx1(line)
        #    if line[0:6] == 'SCALE1': pdb_file.add_scale1(line)
        #    if line[0:6] == 'MTRIX1': pdb_file.add_mtrix1(line)
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM': break
    return
def parse_atoms(args):#models,chains,residues,atoms,debug):
    models = args['models']
    chains = args['chains']
    residues = args['residues']
    atoms = args['atoms']
    debug = args['debug']
    if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  " Parsing permitted atoms. The values supplied:\n" + color.BOLD + "Models:" + color.END + str(models)  + color.BOLD + "\tChains:" + color.END + str(chains) + color.BOLD + "\tResidues:" + color.END + str(residues) + color.BOLD + "\tAtoms:" + color.END + str(atoms)
    if not models:
        if args['macro_ends'] or args['macro_close'] or args['xyz'] or args['traj']: permitted_models = ['*']
        else: permitted_models = [1]
    if models == ['all']: permitted_models = ['*']
    if models and models != ['all']:
        permitted_models = []
        tmp = []
        for part in models: tmp = tmp + part.split(',')
        for part in tmp:
            if len(part.split('-')[0]) == 0: continue
            if len(part.split('-')) == 1: permitted_models.append(int(part))
            if len(part.split('-')) == 2: 
                if int(part.split('-')[0]) < int(part.split('-')[1]):
                    a = int(part.split('-')[0])
                    b = int(part.split('-')[1])
                else:
                    b = int(part.split('-')[0])
                    a = int(part.split('-')[1])
                for k in range(a,b+1): permitted_models.append(k)
    permitted_models = sorted(list(set(permitted_models)))
    if not chains: permitted_chains = ['*']
    else:
        permitted_chains = []
        tmp = []
        for part in chains: tmp = tmp + part.split(',')
        for part in tmp: 
            if part: permitted_chains.append(part)
    permitted_chains = list(set(permitted_chains))
    if not residues: permitted_residues = ['*']
    else:
        permitted_residues = []
        tmp = []
        for part in residues: tmp = tmp + part.split(',')
        for part in tmp:
            if len(part.split('-')[0]) == 0: continue
            if len(part.split('-')) == 1: permitted_residues.append(int(part))
            if len(part.split('-')) == 2: 
                if int(part.split('-')[0]) < int(part.split('-')[1]):
                    a = int(part.split('-')[0])
                    b = int(part.split('-')[1])
                else:
                    b = int(part.split('-')[0])
                    a = int(part.split('-')[1])
                for k in range(a,b+1): permitted_residues.append(k)
    permitted_residues = sorted(list(set(permitted_residues)))
    if not atoms: permitted_atoms = ["CA","C3'"]
    else:
        permitted_atoms = []
        tmp = []
        for part in atoms: tmp = tmp + part.split(',')
        for part in tmp: 
            if part == 'backbone' : permitted_atoms = permitted_atoms + ["CA","C","N","C3'","C4'","C5'","O5'","P","O3'"] 
            elif part == 'all' : 
                permitted_atoms = ['*']
                break
            elif part: permitted_atoms.append(part)
    permitted_atoms = list(set(permitted_atoms))
    if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  " Permitted atoms. The values supplied:\n" + color.BOLD + "Models:" + color.END + str(permitted_models)  + color.BOLD + "\tChains:" + color.END + str(permitted_chains) + color.BOLD + "\tResidues:" + color.END + str(permitted_residues) + color.BOLD + "\tAtoms:" + color.END + str(permitted_atoms)
    return permitted_models, permitted_chains, permitted_residues, permitted_atoms
def gen_chain_name(i):
    if not type(i) is int: return '?'
    if i <= 0 or i > 63: return '?'
    if i <= 26: return chr(i+64)
    if i <= 52: return chr(i+70)
    if i <= 62: return chr(i-5)
    return '?'
def snake():
    try: import pygame, random
    except: 
        print "Sorry, can't do this for you..."
        sys.exit()
    try: import pygame.locals 
    except:
        print "Sorry, can't do this for you..."
        sys.exit()
    def collide(x1, x2, y1, y2, w1, w2, h1, h2):
        if x1+w1>x2 and x1<x2+w2 and y1+h1>y2 and y1<y2+h2:return True
        else:return False
    def die(screen, score):
        f=pygame.font.SysFont('Arial', 30)
        t=f.render('Your score was: '+str(score), True, (0, 0, 0))
        screen.blit(t, (10, 270))
        pygame.display.update()
        pygame.time.wait(2000)
        sys.exit(0)
    def head(x,y,dirs):
        if dirs == 2:
            for k in range(-1,4): pygame.draw.rect(s, (0,0,0), (x+10,y+k*10,9,9)) 
            pygame.draw.rect(s, (0,0,0), (x,y-10,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+20,9,9))
            pygame.draw.rect(s, (0,0,0), (x-10,y+10,9,9))
        if dirs == 3:
            for k in range(-1,4): pygame.draw.rect(s, (0,0,0), (x+k*10,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x-10,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x+10,y-10,9,9))
        if dirs == 0:
            for k in range(-1,4): pygame.draw.rect(s, (0,0,0), (x+10,y-k*10,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y-20,9,9))
            pygame.draw.rect(s, (0,0,0), (x-10,y-10,9,9))
        if dirs == 1:
            for k in range(-1,4): pygame.draw.rect(s, (0,0,0), (x-k*10,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+10,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x-20,y,9,9))
            pygame.draw.rect(s, (0,0,0), (x-10,y-10,9,9))
    def bone(x,y,dirs):
        if dirs == 0:
            for k in range(0,3): pygame.draw.rect(s, (0,0,0), (x,y+k*10,9,9))
            for k in range(1,4): pygame.draw.rect(s, (0,0,0), (x+10,y+k*10,9,9))
        if dirs == 1:
            for k in range(0,3): pygame.draw.rect(s, (0,0,0), (x+k*10,y,9,9))
            for k in range(1,4): pygame.draw.rect(s, (0,0,0), (x+k*10,y+10,9,9))
        if dirs == 2:
            for k in range(0,3): pygame.draw.rect(s, (0,0,0), (x,y-k*10,9,9))
            for k in range(1,4): pygame.draw.rect(s, (0,0,0), (x-10,y-k*10,9,9))
        if dirs == 3:
            for k in range(0,3): pygame.draw.rect(s, (0,0,0), (x-k*10,y,9,9))
            for k in range(1,4): pygame.draw.rect(s, (0,0,0), (x-k*10,y+10,9,9))
    def tail(x,y,dirs):
        if dirs == 0:
            for k in range(0,4): pygame.draw.rect(s, (0,0,0), (x,y+k*10,9,9))
            for k in range(1,6): pygame.draw.rect(s, (0,0,0), (x+10,y+k*10,9,9))
        if dirs == 1:
            for k in range(0,4): pygame.draw.rect(s, (0,0,0), (x+k*10,y,9,9))
            for k in range(1,6): pygame.draw.rect(s, (0,0,0), (x+k*10,y+10,9,9))
        if dirs == 2:
            for k in range(0,4): pygame.draw.rect(s, (0,0,0), (x,y-k*10,9,9))
            for k in range(1,6): pygame.draw.rect(s, (0,0,0), (x-10,y-k*10,9,9))
        if dirs == 3:
            for k in range(0,4): pygame.draw.rect(s, (0,0,0), (x-k*10,y,9,9))
            for k in range(1,6): pygame.draw.rect(s, (0,0,0), (x-k*10,y+10,9,9))
    def number(x,y,k):
        if k == 0:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
                pygame.draw.rect(s, (0,0,0), (x,y+(l+1)*10,9,9))
                pygame.draw.rect(s, (0,0,0), (x+20,y+(l+1)*10,9,9))
        if k == 1:
            for l in range(0,5): pygame.draw.rect(s, (0,0,0), (x+20,y+l*10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+10,y+10,9,9))
        if k == 2:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+20,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+30,9,9))
        if k == 3:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+20,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+30,9,9))
        if k == 4:
            for l in range(0,3): pygame.draw.rect(s, (0,0,0), (x,y+l*10,9,9))
            for l in range(0,5): pygame.draw.rect(s, (0,0,0), (x+20,y+l*10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+10,y+20,9,9))
        if k == 5:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+20,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+30,9,9))
        if k == 6:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+20,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+30,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+30,9,9))
        if k == 7:
            for l in range(0,5): pygame.draw.rect(s, (0,0,0), (x+20,y+l*10,9,9))
            for l in range(0,2): pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
        if k == 8:
            for l in range(0,5): 
                pygame.draw.rect(s, (0,0,0), (x,y+l*10,9,9))
                pygame.draw.rect(s, (0,0,0), (x+20,y+l*10,9,9))
            for l in range(0,3): pygame.draw.rect(s, (0,0,0), (x+10,y+l*20,9,9))
        if k == 9:
            for l in range(0,3): 
                pygame.draw.rect(s, (0,0,0), (x+l*10,y,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+20,9,9))
                pygame.draw.rect(s, (0,0,0), (x+l*10,y+40,9,9))
            pygame.draw.rect(s, (0,0,0), (x,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+10,9,9))
            pygame.draw.rect(s, (0,0,0), (x+20,y+30,9,9))
    def show_score(score,x,y):
        k = 0
        for digit in str("%04d" % score):
            number(x+40*k,y,int(digit))
            k += 1
        return
    def apple(x):
        pygame.draw.rect(s, (0,0,0), (x[0]-10,x[1],9,9))
        pygame.draw.rect(s, (0,0,0), (x[0]+10,x[1],9,9))
        pygame.draw.rect(s, (0,0,0), (x[0],x[1]-10,9,9))
        pygame.draw.rect(s, (0,0,0), (x[0],x[1]+10,9,9))
    xs = [240, 200, 160, 120, 80]
    ys = [165,165,165,165,165]      #23 kolumny, 13 rzedow
    dirs = 1
    score = 0
    pygame.init()
    s=pygame.display.set_mode((850, 500))
    pygame.display.set_caption('Ssssnake...')
    applepos = (15+20*random.randint(0,40), 95+20*random.randint(0,19))
    f = pygame.font.SysFont('Arial', 20)
    clock = pygame.time.Clock()
    while True:
        clock.tick(10)
        for e in pygame.event.get():
            if e.type == pygame.QUIT:
                sys.exit(0)
            elif e.type == pygame.KEYDOWN:
                if e.key == pygame.K_UP and dirs != 0: dirs = 2
                elif e.key == pygame.K_DOWN and dirs != 2: dirs = 0
                elif e.key == pygame.K_LEFT and dirs != 1: dirs = 3
                elif e.key == pygame.K_RIGHT and dirs != 3: dirs = 1
                elif e.key == pygame.K_ESCAPE: sys.exit(0)
        i = len(xs)-1
        while i >= 2:
            if collide(xs[0], xs[i], ys[0], ys[i], 20, 20, 20, 20): die(s, score)
            i-= 1
        if collide(xs[0],applepos[0],ys[0],applepos[1], 20, 20, 20, 20):            
            score+=1;
            xs.append(1000)
            ys.append(1000)
            applepos = (15+20*random.randint(0,40), 95+20*random.randint(0,19))
        if xs[0] < 10: xs[0] = 810 + xs[0]
        if xs[0] > 830: xs[0] = xs[0] - 800
        if ys[0] < 90: ys[0] = 370 + ys[0]
        if ys[0] > 460: ys[0] = ys[0] - 370
        i = len(xs)-1
        while i >= 1:
            xs[i] = xs[i-1]
            ys[i] = ys[i-1]
            i -= 1
        if dirs == 0: ys[0] += 40
        elif dirs == 1: xs[0] += 40
        elif dirs == 2: ys[0] -= 40
        elif dirs == 3: xs[0] -= 40	
        s.fill((255, 255, 255))	
        for k in range(84):
            pygame.draw.rect(s, (0,0,0), (5+10*k,65,9,9))
            pygame.draw.rect(s, (0,0,0), (5+10*k,85,9,9))
            pygame.draw.rect(s, (0,0,0), (5+10*k,485,9,9))
        for k in range(39):
            pygame.draw.rect(s, (0,0,0), (5,95+10*k,9,9))
            pygame.draw.rect(s, (0,0,0), (835,95+10*k,9,9))
        show_score(score,5,5)
        apple(applepos)
        head(xs[0],ys[0],dirs)
        for k in range(1,len(xs)-1): bone(xs[k],ys[k],dirs)
        tail(xs[-1],ys[-1],dirs)
        pygame.display.update()

def informations(pdb_file,args):           # printing PDB info to screen
    debug = args['debug']
    output = args['output'][0]
    if not output: output = 'screen'
    pdb_file.informations(debug)
    return
def read_coordinates(pdb_file,args):  # reading coordinates from PDB file
    debug = args['debug']
    permitted_models, permitted_chains, permitted_residues, permitted_atoms = parse_atoms(args)         # decide, which models/chains/residues/atoms to store
    with open(args['input_file'], 'r') as input_file:
        if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  " File " + args['input_file'] + " opened. Reading the coordinates."
        ###                     # initial values for model
        current_model = 1                   
        store = False
        model_cleaned = True
        model_created = True
        if current_model in permitted_models or '*' in permitted_models:
            pdb_file.create_model(current_model)
            store = True
            if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.GREEN + " Found model " + str(current_model) + color.END
        ###                     # if we have to generate the chain names
        current_chain = 1
        last_index = -99999
        result = ''
        ###
        for line in input_file:
            if line[0:6] == "MODEL ":
                if not model_cleaned and store: 
                    pdb_file.clean(current_model,args['altloc'],permitted_atoms,permitted_residues)
                    result += str(print_coordinates(pdb_file,current_model,args))
                    model_cleaned = True
                current_model = int(line[10:14].strip())
                if (current_model in permitted_models or '*' in permitted_models):
                    store = True
                    if current_model > 1:
                        pdb_file.create_model(current_model)
                        current_chain = 1
                        if type(debug) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + color.GREEN + " Found model " + str(current_model) + color.END
                else: store = False
                model_created = True
            if current_model > 1 and args['bridges']: break
            if line[0:6] == "ENDMDL" and store:
                if not model_cleaned and store: 
                    pdb_file.clean(current_model,args['altloc'],permitted_atoms,permitted_residues)             #cleaning the data
                    result += str(print_coordinates(pdb_file,current_model,args))
                    model_cleaned = True
                    model_created = False
                if args['extended']: break
            if line[0:6] == "ATOM  ":
                if not model_created: 
                    current_model += 1
                    if current_model in permitted_models or '*' in permitted_models:
                        pdb_file.create_model(current_model)
                        current_chain = 1
                        store = True
                        if type(debug) is list and current_model > 1: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.GREEN + " Found model " + str(current_model) + color.END
                    else: store = False
                    model_created = True
                    if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.GREEN + " Found model " + str(current_model) + color.END
                if line[21] == ' ':                  # if there is no chain ID but the resID drops down
                    if int(line[22:26]) < last_index: current_chain += 1
                last_index = int(line[22:26])
                if store: pdb_file.models['mod' + str(current_model)].add_atom(line,current_chain,permitted_chains,permitted_residues,permitted_atoms)
                if model_cleaned: model_cleaned = False
            if line[0:6] == "HETATM":
                if not model_created: 
                    current_model += 1
                    if current_model in permitted_models or '*' in permitted_models:
                        pdb_file.create_model(current_model)
                        current_chain = 1
                        store = True
                        if type(debug) is list and current_model > 1: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.GREEN + " Found model " + str(current_model) + color.END
                    else: store = False
                    model_created = True
                    if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.GREEN + " Found model " + str(current_model) + color.END
                if line[21] == ' ':                  # if there is no chain ID but the resID drops
                    if int(line[22:26]) < last_index: current_chain += 1
                last_index = int(line[22:26])
                if store: pdb_file.models['mod' + str(current_model)].add_atom(line,current_chain,permitted_chains,permitted_residues,permitted_atoms)
                if model_cleaned: model_cleaned = False
            if line[0:6] == "TER   ": current_chain += 1
            if line[0:6] == "END   ":
                if not model_cleaned and store: 
                    pdb_file.clean(current_model,args['altloc'],permitted_atoms,permitted_residues)             #cleaning the data
                    result += str(print_coordinates(pdb_file,current_model,args))
                    model_cleaned = True
    #    if line[0:6] == "ANISOU": pdb_file.add_dbref(line)
    #    if line[0:6] == "CONECT": pdb_file.add_dbref(line)
    #    if line[0:6] == "MASTER": pdb_file.add_dbref(line)
    if not model_cleaned and store: 
        pdb_file.clean(current_model,args['altloc'],permitted_atoms,permitted_residues)             #cleaning the data
        result += str(print_coordinates(pdb_file,current_model,args))
        model_cleaned = True
    if type(debug) is list : print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  " File " + args['input_file'] + " done with the coordinates."
    if args['output'][0] == 'files' : pdb_file.add_ends()
    return result
def print_coordinates(pdb_file,model,args):       # printing coordinates to xyz files
    debug = args['debug']
    output = args['output'][0]
    if type(debug) is list and ('models' in debug or 'chains' in debug or 'residues' in debug or 'atoms' in debug): pdb_file.models['mod' + str(model)].print_model_info(debug)
    if output == 'files':
        pdb_file.models['mod' +str(model)].print_coordinates(pdb_file,args)
        return
    result_pdb = pdb_file.models['mod' + str(model)].print_pdb()[:-1]
    result_xyz = pdb_file.models['mod' + str(model)].print_xyz(args['fourcolumn'],args['traj'],args['columns'],args['hetatoms'])
    if output == 'screen':
        if not args['oxyz']: print result_pdb
        if not args['opdb']: print result_xyz
        return
    if output == 'pipe':
        if not args['oxyz'] and not args['opdb']:
            return result_pdb, result_xyz
        if not args['oxyz'] and args['opdb']:
            return result_pdb
        if not args['opdb'] and args['oxyz']:
            return result_xyz
    if output == 'pipe_xyz':
        return result_xyz
    if output == 'pipe_pdb':
        return result_pdb
    return
def dist(a,b):
    return sqrt(sum([(float(x)-float(y))**2 for x,y in zip(a,b)]))
def macrolink_rev_termini(a):
    if a == 'N': return 'C'
    if a == 'C': return 'N'
    return
def macrolink_ends_distance(perm,orient,distance_matrix):
    distance = float(0)
    for k in range(len(perm)-1):
        distance += distance_matrix[str(perm[k])+str(orient[k]) + '-' + str(perm[k+1])+macrolink_rev_termini(str(orient[k]))]
    distance += distance_matrix[str(perm[-1])+str(orient[-1]) + '-' + str(perm[0])+macrolink_rev_termini(str(orient[0]))]
    return distance
def macrolink_ends_check_distances(loop_arrangement,chains,distance_matrix):
    for k in range(len(loop_arrangement[0])-1):
        d = distance_matrix[str(loop_arrangement[0][k])+str(loop_arrangement[1][k]) + '-' + str(loop_arrangement[0][k+1])+macrolink_rev_termini(str(loop_arrangement[1][k]))]
        if d > 10:
            print 'WARNING!!! It is possible that given chains do not form a macrocomponent. The distance between connected termini is ' + str(d) + ' between ' + chains[k] + ' and ' + chains[k+1] + ' for the best loop!'
    d = distance_matrix[str(loop_arrangement[0][-1])+str(loop_arrangement[1][-1]) + '-' + str(loop_arrangement[0][0])+macrolink_rev_termini(str(loop_arrangement[1][0]))]
    if d > 10:
        print 'WARNING!!! It is possible that given chains do not form a macrocomponent. The distance between connected termini is ' + str(d) + ' between ' + chains[-1] + ' and ' + chains[0] + ' for the best loop!'
    return
def make_macrolink_ends_presentation(loop_arrangements,termini):
    presentation = ''
    for k in range(len(loop_arrangements[0])):
        shift = 0
        if loop_arrangements[1][k] == 'C': shift = 1
        presentation += termini[2*loop_arrangements[0][k]+shift] + '-c-' + termini[2*loop_arrangements[0][k]+(shift+1)%2] + '-b-'
    return presentation[:-3]
def connected_print_xyz(coordinates,atom_list,debug,output='files'):
    result = ''
    for k in range(len(atom_list)):
        if atom_list[k].split('_')[-1] == 'closure':
            result += '99999 closureX closureY closureZ ' + atom_list[k].split('_')[1] + ' closure\n'
            continue
        if len(coordinates[atom_list[k]]) == 3 and atom_list[k].split('_')[-1] != 'closure':
            result += str(k+1) + ' ' + ' '.join(coordinates[atom_list[k]]) + '\n'
        if len(coordinates[atom_list[k]]) != 3 and atom_list[k].split('_')[-1] != 'closure':          
            result += (str(coordinates[atom_list[k]][3]) + ' ' + ' '.join(coordinates[atom_list[k]][:3]) + ' ' + ' '.join(coordinates[atom_list[k]][4:])).strip() + '\n'
    if output == 'screen': print result
    if output == 'files':
        with open('macrolink_ends_result.xyz','w') as myfile: myfile.write(result)
        return
    if output == 'pipe' : 
        return result
    return
def connected_print_pdb(coordinates,atom_list,debug,output='files'):                           #create completelly artificial PDB file, just to watch it in PyMol/Chimera/VMD/etc
    result = ''
    tot = 1 
    for k in range(len(atom_list)):
        if atom_list[k].split('_')[-1] == 'closure':
            result += "ATOM  %(atom)5s  CA  CLS X%(res_nr)4s    closureXclosureYclosureZ  1.00  1.00          C   \n" % {"atom": tot, "res_nr": tot}
            last = [0,0,0]
            if k > 0 : tot += 1
            continue
        if k == 0 :
            last = coordinates[atom_list[k]][:3]
        else:
#            n = int(dist(last,coordinates[atom_list[k]])/2) + 1		#to have the beads in the distance <= 2, for chimera to join them
            n = 1
            delta = [(float(x)-float(y))/n for x,y in zip(coordinates[atom_list[k]],last)]
            res = ''
            if len(coordinates[atom_list[k]]) == 6: res = coordinates[atom_list[k]][5]
            for l in range(n):
                if l == 0 and res == 'CYS': res_name = 'CYS'
#                if l == 0 and res not in amino_acids: res_name = 'ION'
                else: res_name = 'XXX'
                result += "ATOM  %(atom)5s  CA  %(res_name)3s X%(res_nr)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00          C   \n" % {
                                    "atom": l+tot, "res_name": res_name, "res_nr": l+tot,  "x": '{0:.3f}'.format(float(last[0]) + l*delta[0]),
                                    "y": '{0:.3f}'.format(float(last[1]) + l*delta[1]), "z": '{0:.3f}'.format(float(last[2]) + l*delta[2])}
            tot += n
            last = coordinates[atom_list[k]][:3]
    result += "ATOM  %(atom)5s  CA  %(res_name)3s X%(res_nr)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00          C   \n" % {
                                    "atom": tot, "res_name": res_name, "res_nr": tot,  "x": float(last[0]),
                                    "y": float(last[1]), "z": float(last[2])}
    result += "END"
    if output == 'screen': print result
    if output == 'files':
        with open('macrolink_ends_result.pdb','w') as myfile: myfile.write(result)
        return
    if output == 'pipe' : return result
    return
def macrolink_close_find_close(coordinates,chain1,chain2):
    chain1_coords = []
    chain2_coords = []
    bridges = []
    for key in coordinates:
        if '_'.join(key.split('_')[:2]) == chain1: chain1_coords.append([key,coordinates[key][:3]])
        if '_'.join(key.split('_')[:2]) == chain2: chain2_coords.append([key,coordinates[key][:3]])
    for res1 in chain1_coords:
        for res2 in chain2_coords:
            d = dist(res1[1],res2[1])
            if d < 10: 
                bridges.append([res1[0].split('_')[2],res2[0].split('_')[2],d])
    bridges.sort(key=lambda x:x[2])
    return bridges
def build_macro_close(pairs,bridges,iteration): 
    base_shift = int(iteration/len(pairs))
    minor_shift = iteration % len(pairs)
    arrangement = []
    d = 0
    for k in range(len(pairs)):
        shift = base_shift
        if k < minor_shift : shift += 1
        pair_name =  str(pairs[k][0]) + '-' +  str(pairs[k][1])
        if shift >= len(bridges[pair_name]): return [[], 999999]
        arrangement.append([str(pairs[k][0]) + '-' + str(bridges[pair_name][shift][0]),str(pairs[k][1]) + '-' + str(bridges[pair_name][shift][1])])
        d += bridges[pair_name][shift][2]
    return arrangement, d
def macrolink_check_length(macrolink,len_vector):
    if len(macrolink) != len(len_vector): return [False for k in range(len(len_vector))]
    residues = [[] for pair in macrolink]
    for pair in macrolink:
        for k in range(2):
            residues[int(pair[k].split('-')[0])].append(float(pair[k].split('-')[1]))
    vector = [abs(residues[k][0]-residues[k][1]) > len_vector[k] for k in range(len(len_vector))]
    return vector
def macrolink_close_make_presentation(arrangement,chains):
    presentation = ''
    macrolink = []
    for k in range(1,len(arrangement[0])):
        macrolink.append(arrangement[0][k-1][0])
        if arrangement[0][k][1][0] != macrolink[-1][0]: arrangement[0][k] = list(reversed(arrangement[0][k]))
        macrolink.append(arrangement[0][k][1])
    macrolink.append(arrangement[0][-1][0])
    macrolink.append(arrangement[0][0][1])
    for k in range(0,len(macrolink),2):
        a = int(macrolink[k].split('-')[1])
        b = int(macrolink[k+1].split('-')[1])
        c = int(macrolink[k].split('-')[0])
        presentation += chains[c] + '_' + str(a) + '-c-' + chains[c] + '_' + str(b) + '-b-'
    return presentation[:-3]
def print_shortcut(cycle,output='screen'):
    cycle = cycle.split('-')
    k = 1
    while k < len(cycle)-2:
        if cycle[k] == 'c' and cycle[k+2] == 'c':
            cycle.pop(k+1)
            cycle.pop(k)
            k -= 2
        k += 2
    length = len(cycle)
    result = '-'.join(cycle)
    result = re.sub('-b-',' <-> ',re.sub('-c-',' ... ',result))
    if output == 'screen':
        print result
        return
    if output == 'files':
        return
    if output == 'pipe':
        return result,length
    return
def print_macrolink(pdb_file,args):        # printing components chosen automatically from a given set of chains to files
    tmp = []
    coordinates = {}
    graph = {}
    chains = '' 
    output = args['output'][0]
    if args['macro_ends']:
        for x in args['macro_ends']: chains += x + ','
    if args['macro_close']: 
        for x in args['macro_close']: chains += x + ','
    chains = re.sub(',,',',',chains[:-1])
    tmp = re.split(', |,| ',chains)
    k = 0
    while k < len(tmp):
        if tmp[k] == '' :
            tmp.pop(k)
            k -= 1
        k += 1
    for x in tmp:
        if not(x.split('_')[0].isdigit() and len(x.split('_')) == 2 and len(x.split('_')[1]) == 1 and x.split('_')[1].isalpha()):     #check input 
            print "Wrong input format. Expected 'ModelNumber(int)_Chain(alphanumeric,1char)'. Got " + x + ". Exiting."
            sys.exit()
    if not args['columns']: columns = []
    else: columns = args['columns']    
    if args['fourcolumn']: columns = []         # this 'four' is useless actually, but to preserve consistency...
    args['chains'] = tmp
    args['columns'] = ['model','chain','index','x','y','z'] + columns
    args['output'] = ['pipe_xyz']
    result = read_coordinates(PDB,args)
    for atom in result.splitlines(): coordinates['_'.join(atom.split()[:3])] = atom.split()[3:]
    permutations = [[0] + list(perm) for perm in itertools.permutations(range(1,len(tmp)))]
    if args['macro_ends']:
        distance_matrix = {}     # distance matrix
        termini = []
        orientations = [perm for perm in itertools.product('NC',repeat=len(tmp))]
        for chain in tmp: 
            termini.append(chain + '_' + str(pdb_file.models['mod' + chain.split('_')[0]].chains[chain.split('_')[1]].CEnd_number))
            termini.append(chain + '_' + str(pdb_file.models['mod' + chain.split('_')[0]].chains[chain.split('_')[1]].NEnd_number))
        for k in range(len(termini)):
            for l in range(2*int(k/2)):             # to avoid closing begin and end of the same chain
                n1 = int(k/2)
                n2 = int(l/2)
                c1 = 'N'
                c2 = 'N'
                if k % 2 == 0 : c1 = 'C'
                if l % 2 == 0 : c2 = 'C'
                distance_matrix[str(n1)+c1 + '-' + str(n2)+c2] = dist(coordinates[termini[k]][:3],coordinates[termini[l]][:3])
                distance_matrix[str(n2)+c2 + '-' + str(n1)+c1] = dist(coordinates[termini[k]][:3],coordinates[termini[l]][:3])
        test_arrangements = []          # we build every arrangement of the chains and find that one with the smallest distances between termini
        for perm in permutations:
            for orient in orientations: test_arrangements.append([perm,orient,macrolink_ends_distance(perm,orient,distance_matrix)])
        test_arrangements.sort(key=lambda x:x[2])
        macrolink_ends_check_distances(test_arrangements[0],tmp,distance_matrix)
        presentation = make_macrolink_ends_presentation(test_arrangements[0],termini)
    if args['macro_close']:
        bridges = {}
        len_vector = [float(pdb_file.models['mod' + chain.split('_')[0]].chains[chain.split('_')[1]].CEnd_number - pdb_file.models['mod' + chain.split('_')[0]].chains[chain.split('_')[1]].NEnd_number)/10 for chain in chains.split(',')]             #how many residues from each chain we want to include (in this case it is 10)
        for k in range(len(tmp)):
            for l in range(k):
                bridges[str(k) + '-' + str(l)] = macrolink_close_find_close(coordinates,tmp[k],tmp[l])
        test_arrangements = []
        for perm in permutations:
            pairs = []
            for k in range(len(perm)-1):
                pairs.append(list(reversed(sorted([perm[k],perm[k+1]]))))
            pairs.append(list(reversed(sorted([perm[0],perm[-1]]))))
            iteration = 0
            vector = [False for k in range(len(tmp))]
            while all(value != True for value in vector):
                macrolink,d = build_macro_close(pairs,bridges,iteration)
                iteration += 1
                vector = macrolink_check_length(macrolink,len_vector)
                if d > 10*len(len_vector): break
            test_arrangements.append([macrolink,d])
        test_arrangements.sort(key=lambda x:x[1])
        if not test_arrangements or test_arrangements[0] == [[], 999999]: 
            print "ERROR!!! It seems that there is no clear way to connect the chain or the chains supplied are to far away oto form a component. If you desire to connect these chains, please supply the list of inter-chain connections. I cannot guess them for you..."
            return
        presentation = macrolink_close_make_presentation(test_arrangements[0],tmp)
    shortcut = print_shortcut(presentation,output='pipe')[0]
    atom_list = make_atom_list(shortcut,coordinates)
    chains = re.sub(',','_',re.sub('_','',chains))
    print chains + '\t\t' + shortcut
    result_xyz = connected_print_xyz(coordinates,atom_list,debug=args['debug'],output='pipe')
    result_pdb = connected_print_pdb(coordinates,atom_list,debug=args['debug'],output='pipe')
    if output == 'screen': 
        if not args['opdb']:
            print result_xyz
        if not args['opdb'] and not args['oxyz']:
            print '================================================================'
        if not args['oxyz']:
            print result_pdb
        return
    if output == 'files':
        if not args['opdb']:
            with open('macro_' + chains + '.xyz','w') as output_file:
                output_file.write(result_xyz)
        if not args['oxyz']:
            with open('macro_' + chains + '.pdb','w') as output_file:
                output_file.write(result_pdb)
        return
    if output == 'pipe_xyz':
        return result_xyz
    if output == 'pipe_pdb':
        return result_pdb
    if output == 'pipe':
        return result_xyz, result_pdb
    return
def print_commands(PDB,args):                         # printing commands for Wanda's program
    if not args['commands'] or 'no' in args['commands'] : return
    result = ''
    for bridge in PDB.bridge_list:
        if bridge not in PDB.cross_bridge_list and PDB.bridges[bridge].peptide and PDB.bridges[bridge].size >= min_loop_length:
            if 'bridge_type' in args['commands']:
                result += PDB.bridges[bridge].type + ' '
            if 'ready' in args['commands']:
                result += './' + WandaLassoProgram + ' '
            result += args['input_file'] + '_' + PDB.bridges[bridge].res1_chain + '.xyz ' + str(PDB.bridges[bridge].res1_index) + ' ' + str(PDB.bridges[bridge].res2_index) + '\n'
    result = result[:-1]
    if args['output'][0] == 'pipe':
        return result
    print result
    return
def print_cross_bridges(PDB,args):                    # printing cross-chain bridges to screen to screen
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.CYAN + " Printing cross-chains bridges" + color.END
    result = ''
    for bridge in PDB.bridges:
        if bridge in PDB.cross_bridge_list:
            if 'bridge_type' in args['commands']: result += PDB.bridges[bridge].type + ' '
            result += PDB.bridges[bridge].res1_chain + ' ' + str(PDB.bridges[bridge].res1_index) + ' ' + PDB.bridges[bridge].res2_chain + ' ' + str(PDB.bridges[bridge].res2_index) + '\n'
    result = result[:-1]
    if args['output'] == 'pipe':
        return result
    print result
    return

#added by Ola G
def find_chains_in_cycle(cycle):
    result = ''
    for k in range(len(cycle)):
        if cycle.split('_')[1][0] not in result:
            result += cycle.split('_')[1][0]
    return result

def update_shortcut(shortcut,bridges_types):
    result = []
    for atom in shortcut.split():
        if '_' in atom:
            model,chain,res = atom.split('_')
            result.append(bridges_types[atom] + '-' + model + '_' + chain + res)
        else: result.append(atom)
    return ' '.join(result)

def generate_sphere_point(R=1):
    random = list(np.random.random_sample((2,)))
    theta = 2*np.pi*random[0]
    phi = np.arccos(1 - 2*random[1])
    x = str(round(R*np.sin(phi)*np.cos(theta),3))
    y = str(round(R*np.sin(phi)*np.sin(theta),3))
    z = str(round(R*np.cos(phi),3))
    return [x,y,z]

def generate_closure(res1,res2,closure_type,com=False):
    if not com: com = [str(round((float(a)+float(b))/2,3)) for a,b in zip(res1[1:4],res2[1:4])]
    if closure_type == 'direct': return ' '.join([str(round((float(a)+float(b))/2,3)) for a,b in zip(res1[1:4],res2[1:4])])
    if closure_type == 'one_point': return ' '.join([str(float(a)+float(b)) for a,b in zip(com,generate_sphere_point(50))])
    if closure_type == 'two_points': return ' '.join([str(float(a)+float(b)) for a,b in zip(com,generate_sphere_point(50))]) + '\n' + ' '.join([str(float(a)+float(b)) for a,b in zip(com,generate_sphere_point(25))])
    if closure_type == 'parallel':
        direction = generate_sphere_point(50)
        point1 = [float(b) + float(a) for a,b in zip(res1[1:4],direction)]
        point2 = [float(b) + float(a) for a,b in zip(res2[1:4],direction)]
        normalizer = np.srqt(sum([a**2 for a in point1]))
        point1 = [str(round(a/normalizer,3)) for a in point1]
        normalizer = np.srqt(sum([a**2 for a in point2]))
        point2 = [str(round(a/normalizer,3)) for a in point2]
        return ' '.join([str(float(a)+float(b)) for a,b in zip(com,point1)]) + '\n' + ' '.join([str(float(a)+float(b)) for a,b in zip(com,point1)])
    return '0 0 0'

def save_closure_pdb(name,coordinates,args):
    return

def save_closure_xyz(name,coordinates,args):
    result = ''
    result2 = ''
    found = False
    coordinates = coordinates.splitlines()
    for l in range(len(coordinates)):
        if len(coordinates[l].split()) == 6 and coordinates[l].split()[5] == 'closure': 
            found = True
            continue
        if not found : result += coordinates[l] + '\n'
        if found: result2 += coordinates[l] + '\n'
    with open(name,'w') as output: 
        output.write(result2)
        output.write(result)
#    n = args['closure_num'][0]
#    if args['closure'][0] == 'direct': n = 1
#    coordinates = coordinates.splitlines()
#    for k in range(n):
#        result = ''
#        for l in range(len(coordinates)):
#            if len(coordinates[l].split()) == 6 and coordinates[l].split()[5] == 'closure': 
#                res1 = coordinates[l-1].split()
#                res2 = coordinates[(l+1)%len(coordinates)].split()
#                points = generate_closure(res1,res2,args['closure'][0]).splitlines()
#                i = 0
#                for point in points:
#                    result += str(9999-i) + ' ' + point + '\n'
#                    i += 1
#            else: result += coordinates[l] + '\n'
#        with open(name + '_' + str(k) + '.xyz','w') as output: output.write(result)
    return

def print_extended_bridges(pdb_file,args):     
    graph = {}
    coordinates = {}
    chain_residues = {}
    chains = []
    bridges = []
    bridges_residues = {}
    for bridge in pdb_file.bridges:
        if args['noions'] and pdb_file.bridges[bridge].type != 'SS': continue
        chain1 = pdb_file.bridges[bridge].res1_chain
        chain2 = pdb_file.bridges[bridge].res2_chain
    	if args['chains'] != '' and chain1 not in args['chains']: continue
    	if args['chains'] != '' and chain2 not in args['chains']: continue
        res1 = pdb_file.bridges[bridge].res1_index
        res2 = pdb_file.bridges[bridge].res2_index
        if pdb_file.bridges[bridge].res1_type == 'HOH' or pdb_file.bridges[bridge].res2_type == 'HOH': continue
        bridges_residues['1_'+chain1+'_'+str(res1)] = pdb_file.bridges[bridge].res1_type
        bridges_residues['1_'+chain2+'_'+str(res2)] = pdb_file.bridges[bridge].res2_type
        if chain1 == chain2 and abs(res2-res1) == 1: continue
        chains.append('1_'+chain1)
        chains.append('1_'+chain2)
        bridges.append(['1_'+chain1,str(res1),'1_'+chain2,str(res2)])
    chains = list(set(chains))
    cycles = []
    output = args['output'][0]
    if chains:
        args['chains'] = chains
        args['output'] = ['pipe_xyz']
        args['hetatoms'] = True
        if not args['columns']: args['columns'] = ['index']
        args['columns'] = ['model', 'chain', 'index', 'x', 'y', 'z'] + args['columns']
        result = read_coordinates(pdb_file,args)
        for atom in result.splitlines():  coordinates['_'.join(atom.split()[:3])] = atom.split()[3:]
        for bridge in bridges:
             if bridge[0]+'_'+bridge[1] not in coordinates and 'pipe' not in output:
                print "ERROR!!! It seems that there is no residue " + bridges_residues[bridge[0] + '_' + bridge[1]] + '-' + bridge[0] + bridge[1] + ". Please fix this."
                sys.exit()
             if bridge[2]+'_'+bridge[3] not in coordinates and 'pipe' not in output:
                print "ERROR!!! It seems that there is no residue " + bridges_residues[bridge[2] + '_' + bridge[3]] + '-' + bridge[2] + bridge[3] + ". Please fix this."
                sys.exit()
        if args['closure']:                                                     # adding 'bridge' via closure. What type of closure will be decided later.
            for chain in chains: 
                res1 = pdb_file.find_residue(0,chain.split('_')[1],pdb_file.find_chain(0,chain.split('_')[1]).NEnd_number)
                res2 = pdb_file.find_residue(0,chain.split('_')[1],pdb_file.find_chain(0,chain.split('_')[1]).CEnd_number)
                bridges.append([chain,str(res1.name),chain,'closure'])
                bridges_residues[chain+'_'+str(res1.name)] = res1.type
                bridges.append([chain,str(res2.name),chain,'closure'])
                bridges_residues[chain+'_'+str(res2.name)] = res2.type
                bridges_residues[chain+'_closure'] = 'CLS'
        cycles = find_cycles(pdb_file,bridges,coordinates,args,clean='any')
    if len(cycles) == 0 and 'pipe' not in output: 
        print 'ERROR!!! There is no extended loop in this protein.'
        sys.exit()
    atom_list = []
    shortcut = []
    result_xyz = []
    result_pdb = []
    extended_loops = []
    extended_length = []
    code = re.sub('\.pdb','',args['input_file'])
    for k in range(len(cycles)):
        shortcut.append(print_shortcut(cycles[k],output='pipe'))
        atom_list.append(make_atom_list(shortcut[-1][0],coordinates))
        extended_loops.append([k,shortcut[-1][0]])
        extended_length.append([k,shortcut[-1][1]])
        result_xyz.append(connected_print_xyz(coordinates,atom_list[-1],debug=args['debug'],output='pipe'))
        result_pdb.append(connected_print_pdb(coordinates,atom_list[-1],debug=args['debug'],output='pipe'))
    extended_length = sorted(extended_length,key=lambda x: x[1])
    extended_loops = [[extended_length.index([k,l]),extended_loops[k][1]] for k,l in extended_length]
    shortcut = [update_shortcut(shortcut[k][0],bridges_residues) for k,l in extended_length]
    result_xyz = [result_xyz[k] for k,l in extended_length]
    result_pdb = [result_pdb[k] for k,l in extended_length]
    for k in range(len(extended_loops)):
        chh = find_chains_in_cycle(cycles[extended_length[k][0]])
        if 'pipe' not in output: print str(k) + '\t' + code + '.pdb_' + chh + '_' + str(k) + '.xyz' + '\t' + shortcut[k]
        if output == 'screen': 
            if not args['opdb']:
                print result_xyz[k]
            if not args['opdb'] and not args['oxyz']:
                print '================================================================'
            if not args['oxyz']:
                print result_pdb[k]
        if output == 'files':
            if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  ' Saving extended loops to files "' + code + '.pdb_' + chh + '_' + str(k) + '.xyz". Cycles in shorted version are (<-> denotes bridge):'
            if not args['opdb']:
                if 'closure' in shortcut[k]: save_closure_xyz(code + '.pdb_' + chh + '_' + str(k) + '.xyz',result_xyz[k],args)
                else:
                    with open(code + '.pdb_' + chh + '_' + str(k) + '.xyz','w') as output_file:
                        output_file.write(result_xyz[k])
            if not args['opdb']:
                if 'closure' in shortcut[k]: save_closure_pdb(code + '.pdb_' + chh + '_' + str(k) + '.pdb',result_pdb[k],args)
                else:
                    with open(code + '.pdb_' + chh + '_' + str(k) + '.pdb','w') as output_file:
                        output_file.write(result_pdb[k])
    if output == 'pipe_xyz':
        for k in range(len(result_xyz)):
            result_xyz[k] = [data.split() for data in result_xyz[k].splitlines()]
        return extended_loops, result_xyz
    if output == 'pipe_pdb':
        for k in range(len(result_pdb)):
            result_pdb[k] = result_pdb[k].splitlines()
        return extended_loops, result_pdb
    if output == 'pipe':
        for k in range(len(result_xyz)):
            result_xyz[k] = [data.split() for data in result_xyz[k].splitlines()]
        for k in range(len(result_pdb)):
            result_pdb[k] = result_pdb[k].splitlines()
        return extended_loops, result_xyz, result_pdb
    return
def suggest_bridges(PDB,args):                        # printing close deterministic loops for further link analysis to screen
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + color.CYAN + " Suggesting bridges" + color.END
    bridges = []
    possible = []
    for bridge in PDB.bridges:
        if bridge not in PDB.cross_bridge_list and PDB.bridges[bridge].type != 'OTHER': bridges.append(bridge)
    for perm in itertools.combinations(range(len(bridges)),2): possible.append([bridges[perm[0]],bridges[perm[1]]])
    for perm in itertools.combinations(range(len(bridges)),3): possible.append([bridges[perm[0]],bridges[perm[1]],bridges[perm[2]]])
    for perm in itertools.combinations(range(len(bridges)),4): possible.append([bridges[perm[0]],bridges[perm[1]],bridges[perm[2]],bridges[perm[3]]])
    if len(possible) == 0 : return
    result = ''
    name = re.sub('.pdb', '', args['input_file'])
    for bridge in possible:
        order = []
        for k in range(len(bridge)):
            bridgeClass = PDB.bridges[bridge[k]]
            order.append((bridgeClass.res1_index,k))
        order = [y for x,y in sorted(order)]
        for k in order:
            bridgeClass = PDB.bridges[bridge[k]]
            result += name + bridgeClass.res1_chain + ' ' + str(bridgeClass.res1_index) + ' ' + str(bridgeClass.res2_index) + '#'
        result = result[:-1] + '\n'
    if args['output'] == 'pipe':
        return result[:-1]
    print result[:-1]
    return
def close_chains(pdb_file,k,j,args):
    chain1 = pdb_file.models[pdb_file.model_list[0]].chains[pdb_file.models[pdb_file.model_list[0]].chain_list[k]]
    chain2 = pdb_file.models[pdb_file.model_list[0]].chains[pdb_file.models[pdb_file.model_list[0]].chain_list[j]]
    return bool(float(args['factor'][0])*(chain1.R + chain2.R) > sqrt(sum([(x-y)**2 for x, y in zip(chain1.com,chain2.com)])))
def suggest_chains(PDB,args):                         # printing close chains for further link analysis to screen
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.CYAN + " Suggesting chains" + color.END
    chains = []
    for k in range(len(PDB.models[PDB.model_list[0]].chains)):
        for j in range(k):
            if close_chains(PDB,k,j,args): chains.append([k,j])
    k = 0
    while k < len(chains) and len(chains[k])<3:
        l = k + 1
        while l < len(chains):
            if len(chains[k]) == 2 and len(chains[l]) == 2 and len(list(set(chains[k] + chains[l]))) == 3 : chains.append(list(set(chains[k] + chains[l])))
            if sorted([len(chains[k]),len(chains[l])]) == [2,3] and len(list(set(chains[k] + chains[l]))) == 4 : chains.append(list(set(chains[k] + chains[l])))
            l += 1
        k += 1
    if len(chains) == 0 : return
    result = ''
    name = re.sub('.pdb', '', args['input_file'])
    for chain in chains:
        for k in range(len(chain)):
            chainClass = PDB.models[PDB.model_list[0]].chains[PDB.models[PDB.model_list[0]].chain_list[chain[k]]]
            result += name + chainClass.name + ' ' + str(chainClass.NEnd) + ' ' + str(chainClass.CEnd) + '#'
        result = result[:-1] + '\n'
    if args['output'] == 'pipe':
        return result[:-1]
    print result[:-1]
    return
def make_atom_list(cycle,coordinates):
    atom_list = []
    cycle_split = re.sub(' \.\.\. ','|c|',re.sub(' <-> ','|b|',cycle)).split('|')
    if cycle_split[0] in coordinates: atom_list.append(cycle_split[0])
    for k in range(1,len(cycle_split)-1,2):
        if cycle_split[k] == 'c':
            last_res = int(atom_list[-1].split('_')[-1])
            new_res = int(cycle_split[k+1].split('_')[-1])
            chain = '_'.join(cycle_split[k+1].split('_')[:2])
            if last_res < new_res:
                for l in range(last_res+1,new_res+1): 
                    if chain+'_'+str(l) in coordinates: atom_list.append(chain+'_'+str(l))
            if new_res < last_res:
                for l in range(last_res-1,new_res-1,-1): 
                    if chain+'_'+str(l) in coordinates: atom_list.append(chain+'_'+str(l))
        if cycle_split[k] == 'b' : atom_list.append(cycle_split[k+1])
    return atom_list[:-1]
def clean_cycles(cycles,bridges,clean):
    k = 0
    while k < len(cycles):
        if '-b-' not in cycles[k] or '-c-' not in cycles[k]:
            cycles.pop(k)
            k -= 1
        k += 1
    if clean == 'all':
        for bridge in bridges:              # the cycle must have all or any bridges
            k = 0
            while k < len(cycles):
                option1 = bridge[0]+'_'+bridge[1]+'-b-'+bridge[2]+'_'+bridge[3]
                option2 = bridge[2]+'_'+bridge[3]+'-b-'+bridge[0]+'_'+bridge[1]
                if option1 not in cycles[k] and option2 not in cycles[k]:
                    cycles.pop(k)
                    k -= 1
                k += 1
    k = 0
    while k < len(cycles):          # find representant of class of cycles where we collapse cyclic permutations
        l = k+1
        tmp = cycles[k].split('-')[:-1]
        n = len(tmp)
        cyclic = [[tmp[i - j] for i in range(n)] for j in range(n)]
        cyclic = cyclic + [[list(reversed(tmp))[i - j] for i in range(n)] for j in range(n)]
        while l < len(cycles):
            if cycles[l].split('-')[:-1] in cyclic:
                cycles.pop(l)
                l -= 1
            l += 1
        k += 1
    return cycles
def find_cycles(pdb_file,bridges,coordinates,args,clean):
    graph = {}
    chain_residues = {}
    for bridge in bridges:                                     # create bridge part of the graph and check if residues exist
        if bridge[0]+'_'+bridge[1] not in coordinates and bridge[1] != 'closure':
            print "ERROR!!! It seems that there is no residue " + bridge[0] + '_' + bridge[1] + ". Please fix this."
            sys.exit()
        if bridge[2]+'_'+bridge[3] not in coordinates and bridge[3] != 'closure':
            print "ERROR!!! It seems that there is no residue " + bridge[2] + '_' + bridge[3] + ". Please fix this."
            sys.exit()
        if bridge[0]+'_'+bridge[1] not in graph: graph[bridge[0]+'_'+bridge[1]] = [['b',bridge[2]+'_'+bridge[3]]]
        else: graph[bridge[0]+'_'+bridge[1]].append(['b',bridge[2]+'_'+bridge[3]])
        if bridge[2]+'_'+bridge[3] not in graph: graph[bridge[2]+'_'+bridge[3]] = [['b',bridge[0]+'_'+bridge[1]]]
        else: graph[bridge[2]+'_'+bridge[3]].append(['b',bridge[0]+'_'+bridge[1]])
        if bridge[1] != 'closure' and bridge[3] != 'closure':
            d = dist(coordinates[bridge[0]+'_'+bridge[1]][:3],coordinates[bridge[2]+'_'+bridge[3]][:3])
            if d > 10: print "WARNING!!! The bridge " + bridge[0]+bridge[1] + ' <-> ' + bridge[2]+bridge[3] + ' seems to be wrong as its length is equal ' + str(d) + '.'
        if bridge[0] in chain_residues and bridge[1] != 'closure': chain_residues[bridge[0]].append(int(bridge[1]))
        if bridge[0] not in chain_residues and bridge[1] != 'closure': chain_residues[bridge[0]] = [int(bridge[1])]
        if bridge[2] in chain_residues and bridge[3] != 'closure': chain_residues[bridge[2]].append(int(bridge[3]))
        if bridge[2] not in chain_residues and bridge[3] != 'closure': chain_residues[bridge[2]] = [int(bridge[3])]
    for key in chain_residues:                              # add the chain part of the graph (residues connected by the backbone)
        residue_list = sorted(chain_residues[key])
        for k in range(len(residue_list)-1):
            Nterm = pdb_file.find_chain(key.split('_')[0],key.split('_')[1]).NEnd_number
            Cterm = pdb_file.find_chain(key.split('_')[0],key.split('_')[1]).CEnd_number
            if residue_list[k] >= Nterm and residue_list[k+1] >= Nterm and residue_list[k] <= Cterm and residue_list[k+1] <= Cterm:
                graph[key+'_'+str(residue_list[k])].append(['c',key+'_'+str(residue_list[k+1])])
                graph[key+'_'+str(residue_list[k+1])].append(['c',key+'_'+str(residue_list[k])])
    if type(args['debug']) is list and ('merge' in args['debug'] or 'extended' in args['debug'] or 'macro' in args['debug']): 
        print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  ' The dictionary formed during search for loops:'
        for key in list(graph):
            print '\t' + str(key) + '\t' + str(graph[key])
        print '\t________________________________'
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Searching for cycles.'
    cycles = [path for node in graph for path in dfs(graph, node, node)]
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Found ' + str(len(cycles)) + ' cycles. Cleaning cycles.'
    cycles = clean_cycles(cycles,bridges,clean)
    if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Cycles cleaned. Found ' + str(len(cycles)) + ' clean cycles.'
    if type(args['debug']) is list and ('merge' in args['debug'] or 'extended' in args['debug']): 
        print 'These are:'
        for k in range(len(cycles)):
            print cycles[k]
    return cycles
def merge_xyz_component(pdb_file,current_component,coordinates,bridges,output,args):
    if type(args['debug']) is list:
        print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Component ' + str(current_component) + '. I found the following bridges:'
        for bridge in bridges: print bridge[0] + '_' + bridge[1] + ' <-> ' + bridge[2] + '_' + bridge[3]
        print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Creating dictionary.'
    cycles = find_cycles(pdb_file,bridges,coordinates,args,clean='all')
    if len(cycles) > 1: 
        print 'ERROR!!! The bridges given in component ' + str(current_component) + ' form ' + str(len(cycles)) + ' valid components! Please specialize!'
        sys.exit()
    if len(cycles) == 0: 
        print 'ERROR!!! The bridges given in component ' + str(current_component) + ' cannot form a component! Please fix them!'
        sys.exit()
    shortcut = print_shortcut(cycles[0],output='pipe')[0]
    atom_list = make_atom_list(shortcut,coordinates)
    print str(current_component) + '\t\t' + shortcut
    result_xyz = connected_print_xyz(coordinates,atom_list,debug=args['debug'],output='pipe')
    result_pdb = connected_print_pdb(coordinates,atom_list,debug=args['debug'],output='pipe')
    if output == 'screen': 
        if not args['opdb']:
            print result_xyz
        if not args['opdb'] and not args['oxyz']:
            print '================================================================'
        if not args['oxyz']:
            print result_pdb
        return
    if output == 'files':
        if type(args['debug']) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) + ' Saving cycles to files "component_' + str(current_component) +'.xyz". Cycles in shorted version are (<-> denotes bridge):'
        if not args['opdb']:
            with open('component_' + str(current_component) + '.xyz','w') as output_file:
                output_file.write(result_xyz)
        if not args['opdb']:
            with open('component_' + str(current_component) + '.pdb','w') as output_file:
                output_file.write(result_pdb)
        return
    if output == 'pipe_xyz':
        return result_xyz
    if output == 'pipe_pdb':
        return result_pdb
    if output == 'pipe':
        return result_xyz, result_pdb
    return
def merge_xyz(pdb_file,args):                         # printing components for macrolink to files
    bridges = {}
    chains = {}
    coordinates = {}
    output = args['output'][0]
    debug = args['debug']
    if type(debug) is list: print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  color.CYAN + " Preparing macrocomponents with the connections supplied by the user" + color.END
    if not os.path.isfile(args['xyz'][0]):
        print time.strftime("%d %b %H:%M:%S", time.gmtime()) +  " No " + args['xyz'][0] + " file. Cannot do merging of the chain for macrocomponent analysis."
        return
    with open(args['xyz'][0], 'r') as f:
        for line in f:
            if line[0] == '#': 
                current_component = line[1:].split()[0]
                bridges[current_component] = []
                chains[current_component] = []
            else:
                if '_' in line.split()[0] and '_' in line.split()[2] : 
                    bridges[current_component].append(line.split())
                    chains[current_component].append(line.split()[0])
                    chains[current_component].append(line.split()[2])
                if '_' in line.split()[0] and '_' not in line.split()[2] : 
                    bridges[current_component].append([line.split()[0], line.split()[1], '1_' + line.split()[2], line.split()[3]])
                    chains[current_component].append(line.split()[0])
                    chains[current_component].append('1_' + line.split()[2])                    
                if '_' not in line.split()[0] and '_' in line.split()[2] : 
                    bridges[current_component].append(['1_' + line.split()[0], line.split()[1], line.split()[2], line.split()[3]])
                    chains[current_component].append('1_' + line.split()[0])
                    chains[current_component].append(line.split()[2])
                if '_' not in line.split()[0] and '_' not in line.split()[2] : 
                    bridges[current_component].append(['1_' + line.split()[0], line.split()[1], '1_' + line.split()[2], line.split()[3]])
                    chains[current_component].append('1_' + line.split()[0])
                    chains[current_component].append('1_' + line.split()[2])
    chains_total = []
    for key in chains: chains_total += list(set(chains[key]))
    args['chains'] = chains_total
    args['output'] = ['pipe_xyz']
    if not args['columns']: args['columns'] = []
    args['columns'] = ['model', 'chain', 'index', 'x', 'y', 'z'] + args['columns']
    result = read_coordinates(PDB,args)
    for atom in result.splitlines(): coordinates['_'.join(atom.split()[:3])] = atom.split()[3:]
    for key in bridges: merge_xyz_component(pdb_file,key,coordinates,bridges[key],output,args)
    return
def print_summary(pdb_file,output='screen'):          # print summary to screen
    result = ''
    result += color.BOLD + color.UNDERLINE + "Summary:" + color.END + '\n'
    result += color.BOLD + "Total Number of models:\t\t" + color.END + str(pdb_file.model_number) + '\n'
    result += color.BOLD + "Models stored:\t\t\t" + color.END + str(len(pdb_file.model_list)) + '\n'
    result += color.BOLD + "Total number of chains:\t\t" + color.END + str(pdb_file.total_chain_number) + '\n'
    result += color.BOLD + "Total number of res stored:\t" + color.END + str(pdb_file.total_res_number) + '\n'
    result += color.BOLD + "Number of atoms stored:\t\t" + color.END + str(pdb_file.total_atom_number) + '\n'
    result += color.BOLD + "Number of bridges:\t\t" + color.END + str(len(pdb_file.bridge_list)) + '\n'
    result += color.BOLD + "Number of cross-chain bridges:\t" + color.END + str(pdb_file.cross_bridges_number) + '\n'
    if output == 'screen':
        print result
        return
    if output == 'file':
        with open(PDB.code + '.log','a') as myfile:
            myfile.write(result)
    else: return result

################ CLASSES ################
#### Printing cases for debug option
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   BLINK = '\033[5m'
   INVERTED = '\033[7m'
   END = '\033[0m'

#### PDB_File
class PDB_File:
    def __init__(self, args):
        self.name = args['input_file']
        self.header = ''
        self.date = ''
        self.code = ''
        self.obsolete = False
        self.title = ''
        self.split = False
        self.caveat = ''
        self.compound_list = []
        self.compounds = {}
        self.source = ''
        self.keywords = ''
        self.expmethod = []
        self.current_model = 1
        self.model_number = 1
        self.models = {}
        self.model_list = []
        self.total_chain_number = 0
        self.total_res_number = 0
        self.total_atom_number = 0
        self.bridge_list = []
        self.bridges = {}
        self.bridges_cleaned = False
        self.general_chain_list = []
        self.general_chains = {}
        self.missing_residue_models = []
        self.missing_atom_models = []
        self.zero_occupancy_residue_models = []
        self.zero_occupancy_atom_models = []
        self.cross_bridges_number = 0
        self.cross_bridge_list = []
        self.created_PDB = []
    #### Adding part
    def add_header(self,line):
        self.header = line[10:50]
        self.date = line[50:59]
        self.code = line[62:66]
    def add_obsolete(self,line):
        if line[8:10] == '  ' or line[8:10] == ' 1':
            self.obsolete_date = line[11:20]
            self.obsolete = []
        for code in line[31:75].split():
            self.obsolete.append(code)
    def add_title(self,line):
        if line[8:10] == '  ' or line[8:10] == ' 1':
            self.title = line[10:80].rstrip()
        else : self.title += line[10:80].rstrip()
    def add_split(self,line):
        if line[8:10] == '  ' or line[8:10] == ' 1':
            self.split = []
        for code in line[11:80].split():
            self.split.append(code)
    def add_caveat(self,line):
        if line[8:10] == '  ' or line[8:10] == ' 1':
            self.title = line[19:79].rstrip()
        else : self.title += line[19:79].rstrip()
    def add_compound(self,line):
        if line[10:80].lstrip()[0:6] == "MOL_ID":
            self.current_compound = int(line.split()[-1].rstrip(';'))
            self.compound_list.append('comp' + str(self.current_compound))
            self.compounds[self.compound_list[-1]] = Compound(self.current_compound)
        if line[10:80].lstrip()[0:8] == "MOLECULE": self.compounds[self.compound_list[-1]].add_molecule(line[10:80].lstrip())
        if line[10:80].lstrip()[0:5] == "CHAIN": self.compounds[self.compound_list[-1]].add_chains(line[10:80].lstrip())
        if line[10:80].lstrip()[0:8] == "FRAGMENT": self.compounds[self.compound_list[-1]].add_fragment(line[10:80].lstrip())
        if line[10:80].lstrip()[0:7] == "SYNONYM": self.compounds[self.compound_list[-1]].add_synonym(line[10:80].lstrip())
        if line[10:80].lstrip()[0:2] == "EC": self.compounds[self.compound_list[-1]].add_EC(line[10:80].lstrip())
        if line[10:80].lstrip()[0:10] == "ENGINEERED": self.compounds[self.compound_list[-1]].add_engineered(line[10:80].lstrip())
        if line[10:80].lstrip()[0:8] == "MUTATION": self.compounds[self.compound_list[-1]].add_mutation(line[10:80].lstrip())
        if line[10:80].lstrip()[0:13] == "OTHER_DETAILS": self.compounds[self.compound_list[-1]].add_details(line[10:80].lstrip())
    def add_source(self,line):
        if line[10:80].lstrip()[0:6] == "MOL_ID": self.current_compound = int(line.split()[-1].rstrip(';'))
        if line[10:80].lstrip()[0:10] == "SYNTHETIC": self.compounds['comp' + str(self.current_compound)].add_synthetic(line[10:80].lstrip())
        if line[10:80].lstrip()[0:19] == "ORGANISM_SCIENTIFIC": self.compounds['comp' + str(self.current_compound)].add_org_sc(line[10:80].lstrip())
        if line[10:80].lstrip()[0:15] == "ORGANISM_COMMON": self.compounds['comp' + str(self.current_compound)].add_org_com(line[10:80].lstrip())
        if line[10:80].lstrip()[0:14] == "ORGANISM_TAXID": self.compounds['comp' + str(self.current_compound)].add_taxid(line[10:80].lstrip())           #add missing items from source
    def add_keywords(self,line):
        if line[8:10] == '  ' or line[8:10] == ' 1':
            self.keywords = line[10:79].rstrip()
        else: self.keywords += line[10:79].rstrip()
    def add_expmethod(self,line):
        for method in line[10:79].strip().split(';'):
            self.expmethod.append(method)
    def add_nummdl(self,line):
        self.model_number = int(line[10:14].strip())
    def add_bond(self,line):
        number = len(self.bridge_list)
        self.bridges['bridge' + str(number)] = Bridge(number,line)
        self.bridge_list.append('bridge' + str(number))
    def add_missing_residue(self,line):
        if line[10:13] == '   ' and line[15:18].strip() != '' and line[15:18] != 'RES':
            if line[13:19] == 'MODELS' : 
                self.missing_residue_models.append([int(line[19:21].strip()),int(line[22:25].strip())])
                return
            chain = line[19]
            if chain not in self.general_chains: 
                self.general_chain_list.append(chain)
                self.general_chains[chain] = General_chain(chain)
            self.general_chains[chain].add_missing_residue(line)
    def add_missing_atom(self,line):
        if line[10:13] == '   ' and line[15:18].strip() != '' and line[15:18] != 'RES':
            if line[13:19] == 'MODELS' : 
                self.missing_atom_models.append([int(line[19:21].strip()),int(line[22:25].strip())])
                return
            chain = line[19]
            if chain not in self.general_chains: 
                self.general_chain_list.append(chain)
                self.general_chains[chain] = General_chain(chain)
            self.general_chains[chain].add_missing_atom(line)
    def add_zero_occupancy_residue(self,line):
        if line[10:13] == '   ' and line[15:18].strip() != '' and line[15:18] != 'RES':
            if line[13:19] == 'MODELS' : 
                self.zero_occupancy_residue_models.append([int(line[19:21].strip()),int(line[22:25].strip())])
                return
            chain = line[19]
            if chain not in self.general_chains: 
                self.general_chain_list.append(chain)
                self.general_chains[chain] = General_chain(chain)
            self.general_chains[chain].add_zero_occupancy_residue(line)
    def add_zero_occupancy_atom(self,line):
        if line[10:13] == '   ' and line[15:18].strip() != '' and line[15:18] != 'RES':
            if line[13:19] == 'MODELS' : 
                self.zero_occupancy_atom_models.append([int(line[19:21].strip()),int(line[22:25].strip())])
                return
            chain = line[19]
            if chain not in self.general_chains: 
                self.general_chain_list.append(chain)
                self.general_chains[chain] = General_chain(chain)
            self.general_chains[chain].add_zero_occupancy_atom(line)
    def create_model(self,current_model):
        if current_model > self.model_number : self.model_number = current_model
        self.models['mod' + str(current_model)] = Model(current_model,len(self.models))
        self.model_list.append('mod' + str(current_model))
    def add_residue(self,line):
        if line[11] not in self.general_chains: 
            self.general_chain_list.append(line[11])
            self.general_chains[line[11]] = General_chain(line[11])
        self.general_chains[line[11]].add_residue(line)
    def add_ends(self):
        for f in self.created_PDB:
            with open(f,'a') as myfile:
                myfile.write('END')
    #### Searching
    def find_atom(self,model,chain,residue,name):
        print model,chain,residue,name
        if model == 0: return self.models[self.model_list[0]].chains[chain].residues[residue].atoms[name]
        return self.models['mod' + str(model)].chains[chain].residues[residue].atoms[name]
    def find_residue(self,model,chain,residue):
        if model == 0: return self.models[self.model_list[0]].chains[chain].residues[residue]
        return self.models['mod' + str(model)].chains[chain].residues[residue]
    def find_chain(self,model,chain):
        if model == 0: return self.models[self.model_list[0]].chains[chain]
        return self.models['mod' + str(model)].chains[chain]
    def find_missing_residue(self,model,chain,residue_name):
        if not chain in self.general_chain_list: return 'UNK'
        g_chain = self.general_chains[chain]
        # if the missing residue was in remark 465
        for residue in g_chain.residue_list_missing:
            if residue[2] == str(residue_name) and (residue[0] == str(model) or residue[0] == '*'): return residue[1]
        chain = self.find_chain(model,chain)
        # if the residue is in residue list (seqres) and the numbers in seqres correspond to the indices of residues
        residue_shift = chain.residues[chain.residue_list_sorted[0]].name
        if residue_name - residue_shift < len(g_chain.residue_list) and g_chain.residue_list[residue_name - residue_shift - 1] == chain.residues[chain.residue_list_sorted[-1]].type: return g_chain.residue_list[residue_name - residue_shift]
        # if the residue is in residue list (seqres) but the numbers in seqres are shifted
        try: 
            current_residue_list = [chain.residues[chain.residue_list_sorted[k]].type for k in range(residue_name-residue_shift)]
            for k in range(len(g_chain.residue_list) - residue_name + residue_shift):
                residue_list = [g_chain.residue_list[l] for l in range(k,residue_name-residue_shift+k)]
                if current_residue_list == residue_list: 
                    return g_chain.residue_list[residue_name-residue_shift+k]
                    break
        except: return 'UNK'
        # if I cannot tell you, what residue should that be
        return 'UNK'
    #### Analysis
    #### Cleaning part
    def find_total_chain_number(self,current_model):
        self.total_chain_number += self.models['mod' + str(current_model)].total_chain_number
    def find_total_res_number(self,current_model):
        self.total_res_number += self.models['mod' + str(current_model)].total_res_number
    def find_total_atom_number(self,current_model):
        self.total_atom_number += self.models['mod' + str(current_model)].total_atom_number
    def calculate_cross_bridges_number(self):
        number = 0
        for bridge in self.bridge_list:
            if self.bridges[bridge].res1_chain.strip() != self.bridges[bridge].res2_chain.strip() : 
                number += 1
                self.cross_bridge_list.append(bridge)
        self.cross_bridges_number = number
        return 
    def remove_duplicated_bridges(self):
        k = 0
        while k < len(self.bridge_list):
            l = k + 1
            while l < len(self.bridge_list):
                bridge1 = self.bridges[self.bridge_list[k]]
                bridge2 = self.bridges[self.bridge_list[l]]
                if set([bridge1.res1_index,bridge1.res2_index]) == set([bridge2.res1_index,bridge2.res2_index]) and set([bridge1.res1_chain,bridge1.res2_chain]) == set([bridge2.res1_chain,bridge2.res2_chain]):
                    self.bridge_list.pop(l)
                    l -= 1
                l += 1
            k += 1
    def clean(self,current_model,altloc,atoms,residues):
        if not self.bridges_cleaned:
            self.remove_duplicated_bridges()
            self.calculate_cross_bridges_number()
            for l in range(len(self.bridge_list)):
                self.bridges[self.bridge_list[l]].clean(self)
            self.bridges_cleaned = True
        self.models['mod' + str(current_model)].clean(self,altloc,atoms)
        self.find_total_chain_number(current_model)
        self.find_total_res_number(current_model)
        self.find_total_atom_number(current_model)
        return
    #### Printing part
    def informations(self,debug):
        if type(debug) is list and ('pdb' in debug or 'models' in debug or 'chains' in debug or 'residues' in debug or 'atoms' in debug or 'bridges' in debug): 
            print color.BOLD + color.UNDERLINE + "PDB INFORMATION" + color.END
            print color.BOLD + "Header:\t\t" + color.END + self.header
            print color.BOLD + "Date:\t\t" + color.END + self.date
            print color.BOLD + "Code:\t\t" + color.END + self.code
            if not self.obsolete: print color.BOLD + "Obsolete:\t" + color.END + color.GREEN + "No" + color.END
            if self.obsolete: print color.BOLD + "Obsolete:\t" + color.END + color.RED + str(self.obsolete) + color.END
            print color.BOLD + "Title:\t\t" + color.END + self.title
            if not self.split: print color.BOLD + "Split:\t\t" + color.END + color.GREEN + "No" + color.END
            if self.split: print color.BOLD + "Split:\t\t" + color.END + color.RED + str(self.split) + color.END
            if not self.caveat: print color.BOLD + "Caveat:\t\t" + color.END + color.GREEN + "No" + color.END
            if self.caveat: print color.BOLD + "Caveat:\t\t" + color.END + color.RED + str(self.split) + color.END
            print color.BOLD + "Compounds:\t" + color.END + str(len(self.compounds))
            if type(debug) is list and 'compounds' in debug:
                for k in range(len(self.compound_list)):
                    self.compounds[self.compound_list[k]].print_compound_info()
            print color.BOLD + "Keywords:\t" + color.END + self.keywords
            print color.BOLD + "Exp. method:\t" + color.END + str(self.expmethod)
            print color.BOLD + "Models:\t\t" + color.END + str(self.model_number)
            if type(debug) is list and ('models' in debug or 'chains' in debug or 'residues' in debug or 'atoms' in debug):
                for k in range(len(self.model_list)):
                    self.models[self.model_list[k]].print_model_info(debug)
            if type(debug) is list and ('chains' in debug or 'residues' in debug or 'atoms' in debug):
                print color.CYAN + color.BOLD + "\n\t\tGeneral chains:" + color.END
                for k in range(len(self.general_chain_list)):
                    self.general_chains[self.general_chain_list[k]].print_general_chain_info(debug)
            print color.BOLD + "Bridges:\t" + color.END + str(len(self.bridge_list))
            print color.BOLD + "Cross-chain br:\t" + color.END + str(self.cross_bridges_number)
            if type(debug) is list and 'bridges' in debug:
                for k in range(len(self.bridge_list)):
                    self.bridges[self.bridge_list[k]].print_bridge_info()

#### Model
class Model:
    def __init__(self, number, index):
        self.number = number
        self.name = 'mod' + str(number)
        self.index = index
        self.chains = {}
        self.chain_list = []
        self.total_chain_number = 0
        self.total_res_number = 0
        self.total_atom_number = 0
        self.time = number
    #### Adding part
    def add_atom(self,line,current_chain,chains,residues,atoms):            #we create every chain, however only permitted chains are filled with residues (to save space)
        chain_name = line[21]
        if chain_name == ' ': chain_name = gen_chain_name(current_chain)
        if chain_name in chains or '*' in chains or str(self.number)+'_'+chain_name in chains:
            if chain_name not in self.chain_list: 
                self.chain_list.append(chain_name)
                self.chains[chain_name] = Chain(chain_name,self.number,self.index,len(self.chains))
            self.chains[chain_name].add_atom(line,residues,atoms)
    #### Cleaning part
    def clean(self,pdb_file,altloc,atoms):
        for k in range(len(self.chain_list)):
            self.chains[self.chain_list[k]].clean(pdb_file,altloc,atoms)
        self.find_total_chain_number()
        self.find_total_res_number()
        self.find_total_atom_number()
    def find_total_chain_number(self):
        self.total_chain_number = len(self.chain_list)
    def find_total_res_number(self):
        for k in range(len(self.chain_list)):
            self.total_res_number += self.chains[self.chain_list[k]].total_res_number
    def find_total_atom_number(self):
        for k in range(len(self.chain_list)):
            self.total_atom_number += self.chains[self.chain_list[k]].total_atom_number
    #### Printing part
    def print_pdb(self):
        result = ''
        for chain in self.chains:
            result += self.chains[chain].print_pdb()
            result += "\nTER\n"
        if result:
            result = "MODEL    %(model)4s\n" % {"model": self.number} + result
            result += "ENDMDL\n"
        return result
    def print_xyz(self,four,trajectory,column_format,hetatoms):
        result = ''
        for chain in self.chain_list:
            result += self.chains[chain].print_xyz(four,column_format,hetatoms)
            result += "\n"
        if trajectory and result: result = 't ' + str(self.time) + '\n' + result
        return result
    def print_coordinates(self,pdb_file,args):
        if not args['ochains']:
            with open(args['input_file'] + '_extracted.pdb','a') as output_file:
                if args['input_file'] + '_extracted.pdb' not in pdb_file.created_PDB: pdb_file.created_PDB.append(args['input_file'] + '_extracted.pdb')
                output_file.write(self.print_pdb())
        if not args['ogeneral']:
            for chain in self.chains: 
                if not args['oxyz']:
                    model_info = "MODEL    %(model)4s\n" % {"model": self.number}
                    with open(args['input_file'] + '_' + self.chains[chain].name + '.pdb','a') as output_file:
                        if args['input_file'] + '_' + self.chains[chain].name + '.pdb' not in pdb_file.created_PDB: pdb_file.created_PDB.append(args['input_file'] + '_' + self.chains[chain].name + '.pdb')
                        output_file.write(model_info)
                        output_file.write(self.chains[chain].print_pdb())
                        output_file.write('\nENDMDL\n')
                if not args['opdb']:
                    with open(args['input_file'] + '_' + self.chains[chain].name + '.xyz','a') as output_file:
                        if args['traj']: output_file.write('t ' + str(self.time) + '\n')
                        output_file.write(self.chains[chain].print_xyz(args['fourcolumn'],args['columns'],args['hetatoms']))
#                        output_file.write('\n\n')
    def print_model_info(self,debug):
        print color.BOLD + "\tModel:\t\t" + color.END + str(self.number)
        print color.BOLD + "\tTime:\t\t" + color.END + str(self.time)
        print color.BOLD + "\tChains:\t\t" + color.END + str(len(self.chain_list))
        print color.BOLD + "\tChain names:\t" + color.END + str(self.chain_list)
        if type(debug) is list and ('chains' in debug or 'residues' in debug or 'atoms' in debug):
            for k in range(len(self.chain_list)):
                self.chains[self.chain_list[k]].print_chain_info(debug)
        print "\t__________________________________"

#### Chain
class Chain:
    def __init__(self, name, model, model_index, index):
        self.name = name
        self.model = model
        self.model_index = model_index
        self.index = index
        self.residues = {}
        self.residue_list = []
        self.residue_list_sorted = []
        self.residue_missing_list = []
        self.residue_hetatom_list = []
        self.residue_nonpeptide_list = []
        self.total_atom_number = 0
        self.total_res_number = 0
        self.warning = ''
        self.NEnd = 0
        self.CEnd = 0
        self.NEnd_number = 0
        self.CEnd_number = 0
        self.com = [0,0,0]
        self.R = 0
    #### Adding part
    def add_atom(self,line,residues,atoms):
        residue_index = int(line[22:26])
        if residue_index in residues or '*' in residues or line[0:6] == 'HETATM':
            if residue_index not in self.residue_list: 
                self.residue_list.append(residue_index)
                self.residues[residue_index] = Residue(line,self.name,self.model)
            self.residues[residue_index].add_atom(line,atoms)
        return None
    #### Analysis
    def distance(self,res1,res2):
        coord1 = res1.atoms[res1.atom_list[res1.main]].coordinates
        coord2 = res2.atoms[res2.atom_list[res2.main]].coordinates
        return dist(coord1,coord2)
    #### Cleaning part
    def clean(self,pdb_file,altloc,atoms):           
        for k in range(len(self.residue_list)):
            self.residues[self.residue_list[k]].clean(altloc)
        self.sort_residue_list(pdb_file,atoms,altloc)
        self.extract_hetatoms()
        self.extract_nonpeptide()
        self.extract_missing()
        self.extract_termini()
        self.find_total_res_number()
        self.find_total_atom_number()
        self.find_com()
        self.find_R()
        if self.warning: print self.warning  #### todo to poprawic
    def find_com(self):
        com = self.com
        for res in self.residue_list_sorted:
            atom = self.residues[res].atoms[self.residues[res].atom_list[self.residues[res].main]]
            com = [x + float(y) for x, y in zip(com, atom.coordinates)]
        self.com = [x/len(self.residue_list_sorted) for x in com]
    def find_R(self):
        sumR = 0
        for res in self.residue_list_sorted:
            atom = self.residues[res].atoms[self.residues[res].atom_list[self.residues[res].main]]
            sumR += sum([(x-float(y))**2 for x, y in zip(self.com,atom.coordinates)])
        self.R = sqrt(sumR/len(self.residue_list_sorted))
    def extract_hetatoms(self):
        for k in range(len(self.residue_list)):
            if self.residues[self.residue_list[k]].hetatom : self.residue_hetatom_list.append(self.residue_list[k])
    def extract_nonpeptide(self):
        for k in range(len(self.residue_list)):
            if not self.residues[self.residue_list[k]].peptide : self.residue_nonpeptide_list.append(self.residue_list[k])
    def extract_missing(self):
        for k in range(len(self.residue_list)):
            if self.residues[self.residue_list[k]].missing : self.residue_missing_list.append(self.residue_list[k])
    def extract_termini(self):
        for k in range(len(self.residue_list)):
            if self.residues[self.residue_list[k]].NEnd and self.residues[self.residue_list[k]].peptide : self.NEnd = k
            if self.residues[self.residue_list[k]].CEnd and self.residues[self.residue_list[k]].peptide : self.CEnd = k
        self.NEnd_number = int(self.residues[self.residue_list[self.NEnd]].name)
        self.CEnd_number = int(self.residues[self.residue_list[self.CEnd]].name)
    def find_total_res_number(self):
        self.total_res_number = len(self.residue_list)
    def find_total_atom_number(self):
        for k in range(len(self.residue_list)):
            self.total_atom_number += self.residues[self.residue_list[k]].total_atom_number
    def sort_residue_list(self,pdb_file,atoms,altloc):
        k = 0
        while k < len(self.residue_list):                                         #for each residue found in file
            while len(self.residue_list_sorted) == 0:                               #finding the first acceptable residue				
                if len(self.residues[self.residue_list[k]].atom_list):
                    if not self.residues[self.residue_list[k]].hetatom: self.residue_list_sorted.append(self.residue_list[k])       #"atom" residue with CA? Add it.
                    elif self.residues[self.residue_list[k]].type in nucleotides + amino_acids:                                     #"hetatom" residue, but Amino Acid or nucleotide?
                        self.residue_list_sorted.append(self.residue_list[k])                                                   #add it and make it "peptide" = main chain
                        self.residues[self.residue_list_sorted[-1]].peptide = True
                    else:                                                                                                           # so the case which is left is hetatom residue, which is not "normal"
                        for l in range(len(pdb_file.bridge_list)):                                                                       # for example MSE. These should be connected via "LINK" line
                            chain1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_chain
                            chain2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_chain
                            res1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_index
                            res2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_index
                            index1 = int(self.residues[self.residue_list[k]].name)
                            index2 = int(self.residues[self.residue_list[k+1]].name)
                            if (chain1,chain2) == (self.name,self.name) and set((res1,res2)) == set((index1,index2)) and pdb_file.bridges[pdb_file.bridge_list[l]].type != 'OTHER':       # we accept all bridges but "OTHER"
                                self.residue_list_sorted.append(self.residue_list[k])                                                               # if the bridge connect the residue and the next one, 
                                self.residues[self.residue_list_sorted[-1]].peptide = True                                                          # we add it to the list
                                self.residue_hetatom_list.append(self.residue_list[k])
                k += 1
            for l in range(k,len(self.residue_list)):
                if len(self.residues[self.residue_list[l]].atom_list): break                                                                    # in case there is no main atom in res. In such case the resiude is doubled, TODO
            k = l
            if not len(self.residues[self.residue_list[k]].atom_list): break                                                                    # last residues do not have the main atoms
            n = self.residues[self.residue_list[k]].name - self.residues[self.residue_list_sorted[-1]].name                                     # n - sequential distance between residues
            if n == 1 and not self.residues[self.residue_list_sorted[-1]].hetatom and not self.residues[self.residue_list[k]].hetatom:          # n == 1 - no gap.
                d = self.distance(self.residues[self.residue_list_sorted[-1]],self.residues[self.residue_list[k]])
                if (d < 2.0 ): self.warning += 'WARNING!!! Unnatural distance ' + str(d) + ' (< 2.0 A) between residues ' + str(self.residues[self.residue_list_sorted[-1]].name) + ' and ' + str(self.residues[self.residue_list[k]].name) + '.\n'
                if (d > 4.2 ): self.warning += 'WARNING!!! Unnatural distance ' + str(d) + ' (> 4.2 A) between residues ' + str(self.residues[self.residue_list_sorted[-1]].name) + ' and ' + str(self.residues[self.residue_list[k]].name) + '.\n'
                self.residues[self.residue_list_sorted[-1]].CEnd = False
                self.residue_list_sorted.append(self.residue_list[k])
                self.residues[self.residue_list_sorted[-1]].NEnd = False
            if n > 1 and not self.residues[self.residue_list_sorted[-1]].hetatom and not self.residues[self.residue_list[k]].hetatom:           # gap between residues
                if n >= artifact_gap_size: self.warning += 'WARNING!!! In chain ' + self.name + ' there is a gap of length ' + str(n) + ' between residues ' + str(self.residues[self.residue_list_sorted[-1]].name) + ' and ' + str(self.residues[self.residue_list[k]].name) + '.\n'
                res1 = pdb_file.find_residue(self.model,self.name,self.residue_list_sorted[-1])
                res2 = pdb_file.find_residue(self.model,self.name,self.residue_list[k])
                coordinates1 = res1.atoms[res1.atom_list[res1.main]].coordinates
                coordinates2 = res2.atoms[res2.atom_list[res2.main]].coordinates
                delta = [(float(b_i) - float(a_i))/n for a_i, b_i in zip(coordinates1,coordinates2)]
                n0 = int(self.residues[self.residue_list_sorted[-1]].name)
                for j in range(1,n):
                    residue = pdb_file.find_missing_residue(self.model,self.name,n0 + j)
                    line = "ATOM  %(atom)5s  CA  %(resname)3s %(chain)1s%(res_nr)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C  \n" % {
                                "atom": 1, "resname": residue, "chain": self.name, "res_nr": j+n0,
                                "x": '{0:.3f}'.format(float(coordinates1[0])+j*delta[0],3), "y": '{0:.3f}'.format(float(coordinates1[1])+j*delta[1],3),
                                "z": '{0:.3f}'.format(float(coordinates1[2])+j*delta[2],3)}
                    self.residues[n0 + j] = Residue(line,self.name,self.model)
                    self.residues[self.residue_list_sorted[-1]].CEnd = False
                    self.residue_list_sorted.append(n0 + j)
                    self.residues[self.residue_list_sorted[-1]].NEnd = False
                    self.residues[self.residue_list_sorted[-1]].missing = True
                    self.residues[self.residue_list_sorted[-1]].add_atom(line,atoms)
                    self.residues[self.residue_list_sorted[-1]].clean(altloc)
                    self.residue_list.append(n0 + j)
                self.residues[self.residue_list_sorted[-1]].CEnd = False
                self.residue_list_sorted.append(self.residue_list[k])
                self.residues[self.residue_list_sorted[-1]].NEnd = False
            if n == 1 and (self.residues[self.residue_list_sorted[-1]].hetatom or self.residues[self.residue_list[k]].hetatom):
                bridge_exists = False
                for l in range(len(pdb_file.bridge_list)):
                    chain1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_chain
                    chain2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_chain
                    res1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_index
                    res2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_index
                    index1 = int(self.residues[self.residue_list_sorted[-1]].name)
                    index2 = int(self.residues[self.residue_list[k]].name)
                    if (chain1,chain2) == (self.name,self.name) and set((res1,res2)) == set((index1,index2)) and pdb_file.bridges[pdb_file.bridge_list[l]].type != 'OTHER': 
                        bridge_exists = True
                        break
                if bridge_exists:
                    d = self.distance(self.residues[self.residue_list_sorted[-1]],self.residues[self.residue_list[k]])
                    if (d < 2.0 ): self.warning += 'WARNING!!! Unnatural distance ' + str(d) + ' (< 2.0 A) between residues ' + str(self.residues[self.residue_list_sorted[-1]].name) + ' and ' + str(self.residues[self.residue_list[k]].name) + '.\n'
                    if (d > 4.2 ): self.warning += 'WARNING!!! Unnatural distance ' + str(d) + ' (> 4.2 A) between residues ' + str(self.residues[self.residue_list_sorted[-1]].name) + ' and ' + str(self.residues[self.residue_list[k]].name) + '.\n'
                    self.residues[self.residue_list_sorted[-1]].CEnd = False
                    self.residue_list_sorted.append(self.residue_list[k])
                    self.residues[self.residue_list_sorted[-1]].NEnd = False
                    self.residues[self.residue_list_sorted[-1]].peptide = True
            if n > 1 and (self.residues[self.residue_list_sorted[-1]].hetatom and not self.residues[self.residue_list[k]].hetatom):
                bridge_exists = False
                for l in range(len(pdb_file.bridge_list)):
                    chain1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_chain
                    chain2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_chain
                    res1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_index
                    res2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_index
                    index1 = int(self.residues[self.residue_list_sorted[-1]].name)
                    index2 = int(self.residues[self.residue_list[k]].name)
                    if (chain1,chain2) == (self.name,self.name) and set((res1,res2)) == set((index1,index2)) and pdb_file.bridges[pdb_file.bridge_list[l]].type != 'OTHER':
                        bridge_exists = True
                        break
                if bridge_exists:
                    self.residues[self.residue_list_sorted[-1]].CEnd = False
                    self.residue_list_sorted.append(self.residue_list[k])
                    self.residues[self.residue_list_sorted[-1]].NEnd = False
                    self.residues[self.residue_list_sorted[-1]].peptide = True
                else:
                    res1 = pdb_file.find_residue(self.model,self.name,self.residue_list_sorted[-1])
                    res2 = pdb_file.find_residue(self.model,self.name,self.residue_list[k])
                    coordinates1 = res1.atoms[res1.atom_list[res1.main]].coordinates
                    coordinates2 = res2.atoms[res2.atom_list[res2.main]].coordinates
                    delta = [(float(b_i) - float(a_i))/n for a_i, b_i in zip(coordinates1,coordinates2)]
                    n0 = int(self.residues[self.residue_list_sorted[-1]].name)
                    for j in range(1,n):
                        residue = pdb_file.find_missing_residue(self.model,self.name,n0 + j)
                        line = "ATOM  %(atom)5s  CA  %(resname)3s %(chain)1s%(res_nr)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C  \n" % {
                                    "atom": 1, "resname": residue, "chain": self.name, "res_nr": j+n0,
                                    "x": '{0:.3f}'.format(float(coordinates1[0])+j*delta[0],3), "y": '{0:.3f}'.format(float(coordinates1[1])+j*delta[1],3),
                                    "z": '{0:.3f}'.format(float(coordinates1[2])+j*delta[2],3)}
                        self.residues[n0 + j] = Residue(line,self.name,self.model)
                        self.residues[self.residue_list_sorted[-1]].CEnd = False
                        self.residue_list_sorted.append(n0 + j)
                        self.residues[self.residue_list_sorted[-1]].NEnd = False
                        self.residues[self.residue_list_sorted[-1]].missing = True
                        self.residues[self.residue_list_sorted[-1]].add_atom(line,atoms)
                        self.residue_list.append(n0 + j)
                    self.residues[self.residue_list_sorted[-1]].CEnd = False
                    self.residue_list_sorted.append(self.residue_list[k])
                    self.residues[self.residue_list_sorted[-1]].NEnd = False
                    self.residues[self.residue_list_sorted[-1]].peptide = True

            if n > 1 and self.residues[self.residue_list[k]].hetatom:
                bridge_exists = False
                for l in range(len(pdb_file.bridge_list)):
                    chain1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_chain
                    chain2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_chain
                    res1 = pdb_file.bridges[pdb_file.bridge_list[l]].res1_index
                    res2 = pdb_file.bridges[pdb_file.bridge_list[l]].res2_index
                    index1 = int(self.residues[self.residue_list_sorted[-1]].name)
                    index2 = int(self.residues[self.residue_list[k]].name)
                    if (chain1,chain2) == (self.name,self.name) and set((res1,res2)) == set((index1,index2)) and pdb_file.bridges[pdb_file.bridge_list[l]].type != 'OTHER':
                        bridge_exists = True
                        break
                if bridge_exists:
                    self.residues[self.residue_list_sorted[-1]].CEnd = False
                    self.residue_list_sorted.append(self.residue_list[k])
                    self.residues[self.residue_list_sorted[-1]].NEnd = False
                    self.residues[self.residue_list_sorted[-1]].peptide = True
            k += 1
        return 1
    #### Printing part
    def print_pdb(self):
        result = ''
        atom_starting_index = 1
        for k in range(len(self.residue_list_sorted)):
            result += self.residues[self.residue_list_sorted[k]].print_pdb(atom_starting_index,k)
            atom_starting_index += self.residues[self.residue_list_sorted[k]].total_atom_number
        result = result[:-1]
        return result
    def print_xyz(self,four,column_format,hetatoms):
        result = ''
        atom_starting_index = self.NEnd_number
        for res in self.residue_list_sorted:
            result += self.residues[res].print_xyz(atom_starting_index,four,column_format)
            atom_starting_index += self.residues[res].total_atom_number
        if hetatoms:
            for res in self.residue_hetatom_list: 
                atom = self.residues[res].atom_list[self.residues[res].main]
                result += self.residues[res].atoms[atom].print_xyz(atom_starting_index,four,column_format)
        result = result
        return result
    def print_chain_info(self,debug):
        print color.BOLD + "\t\tChain name:\t" + color.END + str(self.name)
        print color.BOLD + "\t\tNEnd:\t\t" + color.END + self.residues[self.residue_list[self.NEnd]].type + '_' + str(self.residues[self.residue_list[self.NEnd]].name)
        print color.BOLD + "\t\tCEnd:\t\t" + color.END + self.residues[self.residue_list[self.CEnd]].type + '_' + str(self.residues[self.residue_list[self.CEnd]].name)
        print color.BOLD + "\t\tResidues num:\t" + color.END + str(len(self.residue_list))
        print color.BOLD + "\t\tPeptide res:\t" + color.END + str(len(self.residue_list_sorted))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug):
            for k in range(len(self.residue_list_sorted)):
                self.residues[self.residue_list_sorted[k]].print_residue_info(debug)
        print color.BOLD + "\t\tNonPeptide res:\t" + color.END + str(len(self.residue_nonpeptide_list))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug):
            for k in range(len(self.residue_nonpeptide_list)):
                self.residues[self.residue_nonpeptide_list[k]].print_residue_info(debug)
        print color.BOLD + "\t\tMissing res:\t" + color.END + str(len(self.residue_missing_list))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug):
            for k in range(len(self.residue_missing_list)):
                self.residues[self.residue_missing_list[k]].print_residue_info(debug)
        print color.BOLD + "\t\tHetAtom res:\t" + color.END + str(len(self.residue_hetatom_list))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug):
            for k in range(len(self.residue_hetatom_list)):
                self.residues[self.residue_hetatom_list[k]].print_residue_info(debug)
        print color.BOLD + "\t\tResidues:" + color.END
        result = '\t\t['
        for k in range(len(self.residue_list_sorted)):
            if self.residues[self.residue_list[k]].missing: result += color.RED + str(self.residues[self.residue_list[k]].name) + color.END + ', '
            elif not self.residues[self.residue_list[k]].peptide: result += color.YELLOW + str(self.residues[self.residue_list[k]].name) + color.END + ', '
            elif self.residues[self.residue_list[k]].hetatom: result += color.PURPLE + str(self.residues[self.residue_list[k]].name) + color.END + ', '
            else: result += str(self.residues[self.residue_list[k]].name) + ', '
        result = result[:-2]
        result += ']'
        print result
        print color.BOLD + "\t\tColor code:\t" + color.END + "Peptide\t\t" + color.RED + "Missing\t\t" + color.END + color.YELLOW + "Nonpeptide\t\t" + color.END + color.PURPLE + "Peptide hetatom" + color.END
        print "\t\t________________________________\n"

#### General chain
class General_chain:
    def __init__(self, name):
        self.name = name
        self.residue_list = []
        self.residue_list_missing = []
        self.zero_occupancy_residues = []
        self.atom_list_missing = []
        self.zero_occupancy_atoms = []
    #### Adding part
    def add_missing_residue(self,line):
        model = line[13:15].strip()
        if model == '': model = '*'
        self.residue_list_missing.append([model,line[15:18],int(line[20:26].strip()),line[26]])
    def add_missing_atom(self,line):
        model = line[13:15].strip()
        if model == '': model = '*'
        atoms = [atom.strip() for atom in line[28:80].split()]
        res = line[20:25].strip()
        if res:
            if res[-1].isalpha():
                altloc = res[-1]
                res = int(res[:-1])
            else: 
                altloc =  ' '
                res = int(res)
        else:
            res = ''
            altloc = ''
        self.atom_list_missing.append([model,line[15:18],res,altloc,atoms])
    def add_zero_occupancy_residue(self,line):
        model = line[13:15].strip()
        if model == '': model = '*'
        self.zero_occupancy_residues.append([model,line[15:18],int(line[20:24].strip()),line[24]])
    def add_zero_occupancy_atom(self,line):
        model = line[13:15].strip()
        if model == '': model = '*'
        atoms = [atom.strip() for atom in line[28:80].split()]
        res = line[20:25].strip()
        if res[-1].isalpha():
            altloc = res[-1]
            res = int(res[:-1])
        else: 
            altloc =  ' '
            res = int(res)
        self.zero_occupancy_atoms.append([model,line[15:18],res,altloc,atoms])
    def add_residue(self,line):
        for res in line[19:80].split():
            self.residue_list.append(res)
    #### Cleaning part
    #### Printing part
    def print_general_chain_info(self,debug):
        print color.BOLD + "\t\tChain name:\t" + color.END + str(self.name)
        print color.BOLD + "\t\tResidues num:\t" + color.END + str(len(self.residue_list))
        print color.BOLD + "\t\tResidues:\t" + color.END + str(len(self.residue_list))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug): print "\t\t\t" + str(self.residue_list)
        print color.BOLD + "\t\tMissing res:\t" + color.END + str(len(self.residue_list_missing))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug): print "\t\t\t" + str(self.residue_list_missing)
        print color.BOLD + "\t\tZero occ res:\t" + color.END + str(len(self.zero_occupancy_residues))
        if type(debug) is list and ('residues' in debug or 'atoms' in debug): print "\t\t\t" + str(self.zero_occupancy_residues)
        print color.BOLD + "\t\tMissing atoms:\t" + color.END + str(len(self.atom_list_missing))
        if type(debug) is list and ('atoms' in debug): print "\t\t\t" + str(self.atom_list_missing)
        print color.BOLD + "\t\tZero occ atoms:\t" + color.END + str(len(self.zero_occupancy_atoms))
        if type(debug) is list and ('atoms' in debug): print "\t\t\t" + str(self.zero_occupancy_atoms)
        print "\t\t________________________________\n"


#### Residue
class Residue:
    def __init__(self, line, chain, model):
        self.name = int(line[22:26])
        self.chain = chain
        self.model = model
        self.type = line[17:20].strip()
        if self.type in amino_acids + nucleotides: self.standard = True
        else: self.standard = False
        self.atoms = {}
        self.atom_list = []
        self.atom_list_sorted = []
        self.total_atom_number = 0
        self.missing = False
        if line[0:6] == 'ATOM  ': self.hetatom = False
        if line[0:6] == 'HETATM': self.hetatom = True
        self.peptide = not self.hetatom
        self.main = 0
        self.CEnd = True
        self.NEnd = True
    #### Adding part
    def add_atom(self,line,atoms):
        if line[12:16].strip() in atoms or '*' in atoms or line[0:6] == 'HETATM':
            self.atom_list.append(line[12:16].strip() + line[16].strip())
            self.atoms[line[12:16].strip() + line[16].strip()] = Atom(line, self.name, self.chain, self.model)
    #### Cleaning part
    def clean(self,altloc):
        if type(altloc) is list: altloc = altloc[0]
        else: altloc = ''
        self.find_main(altloc)
        self.clean_atom_list(altloc)
        self.find_total_atom_number()
    def find_main(self,altloc):            # get the index of main atom in the list of atoms
        if len(self.atoms) == 0: 
            self.peptide = False
            return
        for k in range(len(self.atom_list)):
            if self.atoms[self.atom_list[k]].name.strip() == 'CA' and self.atoms[self.atom_list[k]].altloc.strip() == altloc: 
                self.main = k
                break
            if self.atoms[self.atom_list[k]].name.strip() == 'CA' and self.atoms[self.atom_list[k]].altloc.strip() == '':
                self.main = k
                break
            if self.atoms[self.atom_list[k]].name.strip() == 'CA': 
                self.main = k
                break
            if self.atoms[self.atom_list[k]].name.strip() == "C3'" and self.atoms[self.atom_list[k]].altloc.strip() == altloc: 
                self.main = k
                break
            if self.atoms[self.atom_list[k]].name.strip() == "C3'" and self.atoms[self.atom_list[k]].altloc.strip() == '':
                self.main = k
                break
            if self.atoms[self.atom_list[k]].name.strip() == "C3'":
                self.main = k
                break
        self.atoms[self.atom_list[self.main]].main = True
    def clean_atom_list(self,altloc):
        for atom in self.atom_list:
            if self.atoms[atom].altloc.strip() == '' or self.atoms[atom].altloc.strip() == altloc or self.atoms[atom].main: 
                if not self.hetatom: self.atom_list_sorted.append(atom)
                if self.hetatom:                                                                                                #for hetatoms we add only the first atom
                    self.atom_list_sorted.append(atom)
                    break
    def find_total_atom_number(self):
        self.total_atom_number = len(self.atom_list_sorted)
    #### Printing part
    def print_pdb(self,atom_starting_index,res_nr):
        result = ''
        for k in range(len(self.atom_list_sorted)):
            result += self.atoms[self.atom_list[k]].print_pdb(atom_starting_index+k,res_nr)
        return result
    def print_xyz(self,atom_starting_index,four,column_format):
        result = ''
        for k in range(len(self.atom_list_sorted)):
            result += self.atoms[self.atom_list_sorted[k]].print_xyz(atom_starting_index+k,four,column_format)
        return result
    def print_point_coordinates(self):
        return self.atoms[self.atom_list[self.main]].coordinates
    def print_residue_info(self,debug):
        print color.BOLD + "\t\t\tResidue name:\t" + color.END + str(self.name)
        print color.BOLD + "\t\t\tResidues type:\t" + color.END + self.type
        if self.standard: print color.BOLD + "\t\t\tStandard:\t" + color.END + color.GREEN + str(self.standard) + color.END
        if not self.standard: print color.BOLD + "\t\t\tStandard:\t" + color.END + color.RED + str(self.standard) + color.END
        print color.BOLD + "\t\t\tNum of atoms:\t" + color.END + str(len(self.atom_list))
        if not self.missing: print color.BOLD + "\t\t\tMissing:\t" + color.END + color.GREEN + str(self.missing) + color.END
        if self.missing: print color.BOLD + "\t\t\tMissing:\t" + color.END + color.RED + str(self.missing) + color.END
        if not self.hetatom: print color.BOLD + "\t\t\tHeteroatom:\t" + color.END + color.GREEN + str(self.hetatom) + color.END
        if self.hetatom: print color.BOLD + "\t\t\tHeteroatom:\t" + color.END + color.RED + str(self.hetatom) + color.END
        if not self.peptide: print color.BOLD + "\t\t\tPeptide:\t" + color.END + color.RED + str(self.peptide) + color.END
        if self.peptide: print color.BOLD + "\t\t\tPeptide:\t" + color.END + color.GREEN + str(self.peptide) + color.END
        if type(debug) is list and 'atoms' in debug:
            for k in range(len(self.atom_list)):
                self.atoms[self.atom_list[k]].print_atom_info(debug)
        print "\t\t\t________________________________\n"


#### Atom
class Atom:
    def __init__(self, line, residue, chain, model):
        self.name = line[12:16]
        self.serial = int(line[6:11])
        self.altloc = line[16].strip()
        self.resname = line[17:20].strip()
        self.insertion = line[26].strip()
        self.coordinates = [line[30:38].strip(),line[38:46].strip(),line[46:54].strip()]
        self.occupancy = line[54:60].strip()
        self.bFactor = line[60:66].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        self.residue = residue
        self.chain = chain
        self.model = model
        self.connect = []
        self.main = False
    #### Cleaning part
    #### Printing part
    def print_pdb(self,atom_index,res_index):
        return "ATOM  %(atom)5s %(type)4s %(resname)3s %(chain)1s%(res_nr)4s    %(x)8s%(y)8s%(z)8s%(occupancy)6s%(bFactor)6s          %(element)2s%(charge)2s\n" % {
                                "atom": atom_index, "type": self.name, "resname": self.resname,
                                "chain": self.chain, "res_nr": self.residue,  "x": self.coordinates[0],
                                "y": self.coordinates[1], "z": self.coordinates[2], "occupancy": self.occupancy,
                                "bFactor": self.bFactor, "element": self.element, "charge": self.charge}
    def print_xyz(self,res_index,four,column_format):
        if column_format: 
            result = ''
            columns = ''
            for column in column_format: columns += column + ' '
            columns = re.split(',| ',columns)
            for column in columns:
                if column == '': continue
                if column == 'index': result += str(self.residue) + " "
                if column == 'x': result += str(self.coordinates[0]) + " "
                if column == 'y': result += str(self.coordinates[1]) + " "
                if column == 'z': result += str(self.coordinates[2]) + " "
                if column == 'chain': result += str(self.chain) + " "
                if column == 'residue': result += str(self.resname) + " "
                if column == 'model': result += str(self.model) + " "
                if column == 'bfactor': result += str(self.bFactor) + " "
                if column == 'occupancy': result += str(self.occupancy) + " "
                if column == 'element': result += str(self.element) + " "
                if column == 'residue_index': result += str(self.residue) + " "
                if column == 'charge': result += str(self.charge) + " "
                if column == 'name': result += str(self.name) + " "
                if column == 'insertion': result += str(self.insertion) + " "
                if column == 'alt_location': result += str(self.altlog) + " "
                if column == '|': result += "| "
            result += '\n'
            return result
        if four: return str(res_index) + " " + str(self.coordinates[0]) + " " + str(self.coordinates[1]) + " " + str(self.coordinates[2]) + "\n"
        else: return str(res_index) + " " + str(self.coordinates[0]) + " " + str(self.coordinates[1]) + " " + str(self.coordinates[2]) + " " + self.resname + "\n"
    def print_atom_info(self,debug):
        print color.BOLD + "\t\t\t\tName:\t\t" + color.END + self.name.strip()
        print color.BOLD + "\t\t\t\tSerial:\t\t" + color.END + str(self.serial)
        print color.BOLD + "\t\t\t\tAlt Location:\t" + color.END + self.altloc
        print color.BOLD + "\t\t\t\tInsertion:\t" + color.END + self.insertion
        print color.BOLD + "\t\t\t\tCoordinates:\t" + color.END + str(self.coordinates)
        print color.BOLD + "\t\t\t\tOccupancy:\t" + color.END + str(self.occupancy)
        print color.BOLD + "\t\t\t\tBeta Factor:\t" + color.END + str(self.bFactor)
        print color.BOLD + "\t\t\t\tElement:\t" + color.END + self.element
        print color.BOLD + "\t\t\t\tCharge:\t" + color.END + self.charge
        if self.main : print color.BOLD + "\t\t\t\tMain:\t\t" + color.END + color.GREEN + str(self.main) + color.END
        if not self.main : print color.BOLD + "\t\t\t\tMain:\t\t" + color.END + color.RED + str(self.main) + color.END
        print "\t\t\t\t________________________________\n"


#### Compound
class Compound:
    def __init__(self, number):
        self.number = number
        self.name = 'comp' + str(number)
        self.molecule = ''
        self.chains = []
        self.fragment = ''
        self.synonym = ''
        self.EC = ''
        self.engineered = False
        self.mutation = False
        self.details = ''
        self.synthetic = False
        self.org_sc = ''
        self.org_com = ''
        self.taxid = ''
    def add_molecule(self,line):
        self.molecule = re.sub('^.*MOLECULE:','',line).strip()
    def add_chains(self,line):
        for char in re.sub('^.*CHAIN:','',line).strip(' ;').split():
            self.chains.append(char.strip(','))
    def add_fragment(self,line):
        self.fragment = re.sub('^.*FRAGMENT:','',line).strip()
    def add_synonym(self,line):
        self.synonym = re.sub('^.*SYNONYM:','',line).strip()
    def add_EC(self,line):
        self.EC = re.sub('^.*EC:','',line).strip()
    def add_engineered(self,line):
        if line.split()[-1].strip(';') == 'YES' : self.engineered = True
    def add_mutation(self,line):
        if line.split()[-1].strip(';') == 'YES' : self.engineered = True
    def add_details(self,line):
        self.details = re.sub('^.*DETAILS:','',line).strip()
    def add_synthetic(self,line):
        if line.split()[-1].strip(';') == 'YES' : self.synthetic = True
    def add_org_sc(self,line):
        self.org_sc = re.sub('^.*ORGANISM_SCIENTIFIC:','',line).strip()
    def add_org_com(self,line):
        self.org_com = re.sub('^.*ORGANISM_COMMON:','',line).strip()
    def add_taxid(self,line):
        self.taxid = re.sub('^.*ORGANISM_TAXID:','',line).strip()
    #### Cleaning part
    #### Printing part
    def print_compound_info(self):
        print color.BOLD + "\tName:\t\t" + color.END + self.name
        print color.BOLD + "\tMolecule:\t" + color.END + self.molecule
        print color.BOLD + "\tChains:\t\t" + color.END + str(self.chains)
        print color.BOLD + "\tFragment:\t" + color.END + self.fragment
        print color.BOLD + "\tSynonym:\t" + color.END + self.synonym
        print color.BOLD + "\tEC:\t\t" + color.END + self.EC
        print color.BOLD + "\tEngineered:\t" + color.END + str(self.engineered)
        print color.BOLD + "\tMutation:\t" + color.END + str(self.mutation)
        print color.BOLD + "\tSynthetic:\t" + color.END + str(self.synthetic)
        print color.BOLD + "\tOrganism Sci:\t" + color.END + self.org_sc
        print color.BOLD + "\tOrganism Com:\t" + color.END + self.org_com
        print color.BOLD + "\tTAXID:\t\t" + color.END + self.taxid
        print color.BOLD + "\tDetails:\t" + color.END + self.details
        print "\t__________________________________"

#### Bridge
class Bridge:
    #### Adding
    def __init__(self, number, line):
        self.number = number
        self.name = 'bridge' + str(number)
        if line[0:6] == 'SSBOND':
            self.serial = int(line[7:10])
            self.res1_index = int(line[17:21])
            self.res1_type = 'CYS'
            self.res1_chain = line[15].strip()
            self.res1_insert = line[21].strip()
            self.res1_atom = 'SG'
            self.res2_index = int(line[31:35])
            self.res2_type = 'CYS'
            self.res2_chain = line[29].strip()
            self.res2_insert = line[35].strip()
            self.res2_atom = 'SG'
            self.type = 'SS'
            self.peptide = True
        if line[0:6] == 'LINK  ':
            self.serial = ''
            self.res1_index = int(line[22:26])
            self.res1_type = line[17:20].strip()
            self.res1_chain = line[21].strip()
            self.res1_insert = line[26].strip()
            self.res1_atom = line[12:16]
            self.res2_index = int(line[52:56])
            self.res2_type = line[47:50].strip()
            self.res2_chain = line[51].strip()
            self.res2_insert = line[56].strip()
            self.res2_atom = line[42:46]
            self.type = ''
            self.peptide = False
        self.res1_NEnd = False
        self.res1_CEnd = False
        self.res2_NEnd = False
        self.res2_CEnd = False
        self.distance = True
        self.res1_symop = line[59:65].strip()
        self.res2_symop = line[66:72].strip()
        try: self.dist = float(line[73:78])
        except ValueError: self.dist = 999
        self.size = abs(self.res1_index - self.res2_index)
        self.exists = True
        self.bond_type()
    #### Cleaning
    def clean(self,pdb_file):
        self.decide_exists(pdb_file)
        if self.exists: self.decide_termini(pdb_file)
        if self.exists: self.decide_peptide(pdb_file)
        self.decide_distance()
    def decide_distance(self):
        if self.dist < 2.0 or self.dist > 4.2: self.distance = False
    def decide_termini(self,pdb_file):
        chain1 = pdb_file.find_chain(0,self.res1_chain)
        chain2 = pdb_file.find_chain(0,self.res2_chain)
        self.res1_CEnd = pdb_file.find_residue(0,self.res1_chain,self.res1_index).CEnd
        self.res2_CEnd = pdb_file.find_residue(0,self.res2_chain,self.res2_index).CEnd
        self.res1_NEnd = pdb_file.find_residue(0,self.res1_chain,self.res1_index).NEnd
        self.res2_NEnd = pdb_file.find_residue(0,self.res2_chain,self.res2_index).NEnd
    def decide_peptide(self,pdb_file):
        self.peptide = pdb_file.find_residue(0,self.res1_chain,self.res1_index).peptide and pdb_file.find_residue(0,self.res2_chain,self.res2_index).peptide
    def decide_exists(self,pdb_file):
        try: res1 = pdb_file.find_residue(0,self.res1_chain,self.res1_index)
        except KeyError: 
            self.exists = False
            return
        try: res2 = pdb_file.find_residue(0,self.res2_chain,self.res2_index)
        except KeyError: 
            self.exists = False
            return
        self.exists = not (res1.missing or res2.missing)
    def bond_type(self):
        if self.res1_atom.strip()[0] == 'S' and self.res2_atom.strip()[0] == 'S': self.type = 'SS'
        elif self.res1_atom.strip()[0] == 'C' and self.res2_atom.strip()[0] == 'N':
            if self.res1_type == 'GLU' and self.res1_atom.strip() == 'CD' and self.res2_type == 'LYS' and self.res2_atom.strip() == 'NZ' : self.type = 'AMIDE'
            elif self.res1_type == 'ASP' and self.res1_atom.strip() == 'CG' and self.res2_type == 'LYS' and self.res2_atom.strip() == 'NZ' : self.type = 'AMIDE'
            elif self.res1_CEnd and self.res2_type == 'LYS' and self.res2_atom.strip() == 'NZ' : self.type = 'AMIDE'
            elif self.res1_type == 'GLU' and self.res1_atom.strip() == 'CD' and self.res2_NEnd : self.type = 'AMIDE'
            elif self.res1_type == 'ASP' and self.res1_atom.strip() == 'CG' and self.res2_NEnd : self.type = 'AMIDE'
            else: self.type = 'AMIDE-like'
        elif self.res1_atom.strip()[0] == 'N' and self.res2_atom.strip()[0] == 'C':
            if self.res1_type == 'LYS' and self.res1_atom.strip() == 'NZ' and self.res2_type == 'GLU' and self.res2_atom.strip() == 'CD' : self.type = 'AMIDE'
            elif self.res1_type == 'LYS' and self.res1_atom.strip() == 'NZ' and self.res2_type == 'ASP' and self.res2_atom.strip() == 'CG' : self.type = 'AMIDE'
            elif self.res1_NEnd and self.res2_type == 'GLU' and self.res2_atom.strip() == 'CD' : self.type = 'AMIDE'
            elif self.res1_NEnd and self.res2_type == 'ASP' and self.res2_atom.strip() == 'CG' : self.type = 'AMIDE'
            elif self.res1_type == 'LYS' and self.res1_atom.strip() == 'NZ' and self.res2_CEnd : self.type = 'AMIDE'
            else: self.type = 'AMIDE-like'
        elif self.res1_atom.strip()[0] == 'C' and self.res2_atom.strip()[0] == 'O':
            if self.res1_type == 'GLU' and self.res1_atom.strip() == 'CD' and self.res2_type == 'THR' and self.res2_atom.strip() == 'OG1' : self.type = 'ESTER'
            elif self.res1_type == 'GLU' and self.res1_atom.strip() == 'CD' and self.res2_type == 'SER' and self.res2_atom.strip() == 'OG' : self.type = 'ESTER'
            elif self.res1_type == 'ASP' and self.res1_atom.strip() == 'CG' and self.res2_type == 'THR' and self.res2_atom.strip() == 'OG1' : self.type = 'ESTER'
            elif self.res1_type == 'ASP' and self.res1_atom.strip() == 'CG' and self.res2_type == 'SER' and self.res2_atom.strip() == 'OG' : self.type = 'ESTER'
            elif self.res1_CEnd and self.res2_type == 'SER' and self.res2_atom.strip() == 'OG' : self.type = 'ESTER'
            elif self.res1_CEnd and self.res2_type == 'THR' and self.res2_atom.strip() == 'OG1' : self.type = 'ESTER'
            else: self.type = 'ESTER-like'
        elif self.res1_atom.strip()[0] == 'O' and self.res2_atom.strip()[0] == 'C':
            if self.res1_type == 'THR' and self.res1_atom.strip() == 'OG1' and self.res2_type == 'GLU' and self.res2_atom.strip() == 'CD' : self.type = 'ESTER'
            elif self.res1_type == 'THR' and self.res1_atom.strip() == 'OG1' and self.res2_type == 'ASP' and self.res2_atom.strip() == 'CG' : self.type = 'ESTER'
            elif self.res1_type == 'SER' and self.res1_atom.strip() == 'OG' and self.res2_type == 'GLU' and self.res2_atom.strip() == 'CD' : self.type = 'ESTER'
            elif self.res1_type == 'SER' and self.res1_atom.strip() == 'OG' and self.res2_type == 'ASP' and self.res2_atom.strip() == 'CG' : self.type = 'ESTER'
            elif self.res1_type == 'THR' and self.res1_atom.strip() == 'OG1' and self.res2_CEnd : self.type = 'ESTER'
            elif self.res1_type == 'SER' and self.res1_atom.strip() == 'OG' and self.res2_CEnd : self.type = 'ESTER'
            else: self.type = 'ESTER-like'
        elif self.res1_atom.strip()[0] == 'C' and self.res2_atom.strip()[0] == 'S':
            if self.res1_type == 'GLU' and self.res1_atom.strip() == 'CD' and self.res2_type == 'CYS' and self.res2_atom.strip() == 'SG' : self.type = 'THIOESTER'
            elif self.res1_type == 'ASP' and self.res1_atom.strip() == 'CG' and self.res2_type == 'CYS' and self.res2_atom.strip() == 'SG' : self.type = 'THIOESTER'
            elif self.res1_CEnd and self.res2_type == 'CYS' and self.res2_atom.strip() == 'SG' : self.type = 'THIOESTER'
            elif self.res1_CEnd and self.res2_type == 'CYS' and self.res2_atom.strip() == 'SG' : self.type = 'THIOESTER'
            else: self.type = 'THIOESTER-like'
        elif self.res1_atom.strip()[0] == 'S' and self.res2_atom.strip()[0] == 'C':
            if self.res1_type == 'CYS' and self.res1_atom.strip() == 'SG' and self.res2_type == 'GLU' and self.res2_atom.strip() == 'CD' : self.type = 'THIOESTER'
            elif self.res1_type == 'CYS' and self.res1_atom.strip() == 'SG' and self.res2_type == 'ASP' and self.res2_atom.strip() == 'CD' : self.type = 'THIOESTER'
            elif self.res1_type == 'CYS' and self.res1_atom.strip() == 'SG' and self.res2_CEnd : self.type = 'THIOESTER'
            elif self.res1_type == 'CYS' and self.res1_atom.strip() == 'SG' and self.res2_CEnd : self.type = 'THIOESTER'
            else: self.type = 'THIOESTER-like'
        else: self.type = 'OTHER'

    #### Printing
    def print_bridge_info(self):
        print color.BOLD + "\tRes1 index:\t" + color.END + str(self.res1_index) + color.BOLD + "\tRes2 index:\t" + color.END + str(self.res2_index)
        print color.BOLD + "\tRes1 type:\t" + color.END + self.res1_type + color.BOLD + "\tRes2 type:\t" + color.END + self.res2_type
        print color.BOLD + "\tRes1 chain:\t" + color.END + self.res1_chain + color.BOLD + "\tRes2 chain:\t" + color.END + self.res2_chain
        print color.BOLD + "\tRes1 ins:\t" + color.END + self.res1_insert + color.BOLD + "\tRes2 ins:\t" + color.END + self.res2_insert
        print color.BOLD + "\tRes1 symmetry:\t" + color.END + self.res1_symop + color.BOLD + "\tRes2 symmetry:\t" + color.END + self.res2_symop
        print color.BOLD + "\tRes1 atom:\t" + color.END + self.res1_atom.strip() + color.BOLD + "\tRes2 atom:\t" + color.END + self.res2_atom.strip()
        print color.BOLD + "\tRes1 NEnd:\t" + color.END + str(self.res1_NEnd) + color.BOLD + "\tRes2 NEnd:\t" + color.END + str(self.res2_NEnd)
        print color.BOLD + "\tRes1 CEnd:\t" + color.END + str(self.res1_CEnd) + color.BOLD + "\tRes2 CEnd:\t" + color.END + str(self.res2_CEnd)
        if self.dist > 4.2 or self.dist < 2.0: print color.BOLD + "\tDistance:\t" + color.END + color.RED + str(self.dist) + color.END
        if self.dist <= 4.2 and self.dist >= 2.0: print color.BOLD + "\tDistance:\t" + color.END + color.GREEN + str(self.dist) + color.END
        if self.type == 'SS': print color.BOLD + "\tType:\t\t" + color.END + color.YELLOW + self.type + color.END
        elif self.type == 'AMIDE' or self.type == 'AMIDE-like': print color.BOLD + "\tType:\t\t" + color.END + color.RED + self.type + color.END
        elif self.type == 'ESTER' or self.type == 'ESTER-like': print color.BOLD + "\tType:\t\t" + color.END + color.GREEN + self.type + color.END
        elif self.type == 'THIOESTER' or self.type == 'THIOESTER-like': print color.BOLD + "\tType:\t\t" + color.END + color.BLUE + self.type + color.END
        else: print color.BOLD + "\tType:\t\t" + color.END + self.type
        print color.BOLD + "\tLoop size:\t" + color.END + str(self.size)
        if self.peptide: print color.BOLD + "\tPeptide:\t" + color.END + color.GREEN + str(self.peptide) + color.END
        else: print color.BOLD + "\tPeptide atoms:\t" + color.END + color.RED + str(self.peptide) + color.END
        if self.exists: print color.BOLD + "\tExists:\t\t" + color.END + color.GREEN + str(self.exists) + color.END
        else: print color.BOLD + "\tAtoms Exist:\t\t" + color.END + color.RED + str(self.exists) + color.END
        print "\t__________________________________"


################ MAIN PART ################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="convert_columns", formatter_class=argparse.RawDescriptionHelpFormatter, description="#################################################################\n\
#   convert_column - script extracting coordinates from PDB     #\n\
#   files and saving them as XYZ file.                          #\n\
#   Script prepares the files for further topological analysis. #\n\
#                                                               #\n\
#   It is module of find_topology script analyzing full         #\n\
#   topology of the protein.                                    #\n\
#                                                               #\n\
#   Date: 02.10.2015, version from " + date + "         	        #\n\
#   Author: " + color.RED + "Pawel Dabrowski-Tumanski" + color.END + "    	                #\n\
#   " + color.RED + "p.dabrowski [at] cent.uw.edu.pl" + color.END + "                             #\n\
#   version 4.0 Beta                                            #\n\
#################################################################")
    parser.add_argument('input_file', action="store", help="The input PDB file") 
    parser.add_argument('-t', '--trajectory', action="store_true", dest="traj", default=False, help="Declare, that the input file is a trajectory")
    parser.add_argument('-f', '--fourcolumn', action="store_true", dest="fourcolumn", default=False, help="Print XYZ output in 4-column format (default 5-column)") 
    parser.add_argument('-b', '--bridges', action="store_true", dest="bridges", default=False, help="Write information about inter-chain bridges.")
    parser.add_argument('-e', '--extended', action="store_true", dest="extended", default=False, help="Write information about extended-type loops (via ions, double bridges etc).")
    parser.add_argument('-s', '--suggest_bridges', action="store_true", dest="sbridge", default=False, help="Suggest bridges for deterministic link analysis based on distance criterium.")
    parser.add_argument('-z', '--suggest_chains', action="store_true", dest="schain", default=False, help="Suggest chains for probabilistic link analysis based on distance criterium.")
    parser.add_argument('-n', '--macro_ends', nargs='*', dest="macro_ends", default=False, help="Prepare xyz files with merged chains. Chains are connected by the termini.")
    parser.add_argument('-m', '--macro_close', nargs='*', dest="macro_close", default=False, help="Prepare xyz files with merged chains. Chains are connected by the closest residues.")
    parser.add_argument('-x', '--xyzfiles', nargs=1, dest="xyz", default=False, help="Prepare xyz files with merged chains. The list of chain connections must be supplied.")
    parser.add_argument('--models', nargs='*', dest="models", default='', help="List of models to be stored.")
    parser.add_argument('--chains', nargs='*', dest="chains", default='', help="List of chains to be stored.")
    parser.add_argument('--residues', nargs='*', dest="residues", default='', help="List of residues to be stored.")
    parser.add_argument('--atoms', nargs='*', dest="atoms", default='', help="List of atom types to be stored.")
    parser.add_argument('--output', nargs=1, dest="output", default=['files'], help="Output of the coordinates.", choices=['files','screen','pipe'])
    parser.add_argument('--columns', nargs='*', dest="columns", default='', help="Format of columns for xyz files. Accepted descriptors are index, x, y, z, name, element, residue, chain, model, residue_index, bfactor, occupancy, charge, insertion, alt_location and '|' as separator. Other values will be omitted.")
    parser.add_argument('--only_general', '--general_only', action="store_true", dest="ogeneral", default=False, help="Print only one file with extracted coordinates chosen.")
    parser.add_argument('--only_chains', '--chains_only', action="store_true", dest="ochains", default=False, help="Print only files for each chain.")
    parser.add_argument('--only_pdb', '--pdb_only', action="store_true", dest="opdb", default=False, help="Print only PDB files.")
    parser.add_argument('--only_xyz', '--xyz_only', action="store_true", dest="oxyz", default=False, help="Print only XYZ files.")
    parser.add_argument('--hetatoms', action="store_true", dest="hetatoms", default=False, help="Print also hetatoms.")
    parser.add_argument('--sug_chain_factor', nargs=1, dest="factor", default=['0.85'], help="Factor used in searching close chains used in link analysis.")
    parser.add_argument('--commands', nargs='*', dest="commands", default=['bridge_type'], choices=['no','bridge_type','ready'], help="Commands format for Wanda lasso program.")
    parser.add_argument('--noions', action="store_true", dest="noions", default=False,  help="If not to include the ions.")
    parser.add_argument('--closure', nargs=1, dest="closure", default=[], help="To close the protein and how.", choices=['one_point','two_points','parallel','direct'])
    parser.add_argument('--closure_num', nargs=1, dest="closure_num", default=100, type=int, help="Number of closures to be done (for two_points, parallel, and one_point).")
    parser.add_argument('--bridge_len', nargs=1, dest="bridge_len", default=[15], help="Maximal bridge length in Angstrems.")
    parser.add_argument('-l', '--altloc', nargs=1, dest="altloc", default=['A'], help="Alternative location of atom to choose. Default 'A'.")
    parser.add_argument('-d', '--debug', nargs='*', choices=['pdb','compounds','models','chains','residues','atoms','coords_pdb','coords_xyz','bridges','extended','merge'], dest="debug", default=False, help=argparse.SUPPRESS)
    parser.add_argument('--version', action='version', version='%(prog)s 4.0')
    args = vars(parser.parse_args())

    PDB = PDB_File(args)
    read_header(PDB,args)                                                                                               # reading header of PDB file
    informations(PDB,args)                                                                                              # printing PDB info to screen
    if not (args['macro_ends'] or args['macro_close'] or args['xyz'] or args['extended']): 
        read_coordinates(PDB,args)                                                                                      # reading coordinates from PDB file and print selected 
        if not args['bridges'] and not args['sbridge'] and not args['schain']: print_commands(PDB,args)                 # printing commands to screen
    if args['macro_ends'] or args['macro_close']: print_macrolink(PDB,args)                                             # printing components chosen automatically from a given set of chains to files
    if args['xyz']: merge_xyz(PDB,args)                                                                                 # printing components for macrolink to files
    if args['extended']: print_extended_bridges(PDB,args)                                                               # printing close chains for further link
    if args['bridges']: print_cross_bridges(PDB,args)                                                                   # printing cross-chain bridges to screen to screen
    if args['sbridge']: suggest_bridges(PDB,args)                                                                       # printing close deterministic loops for further link analysis to screen
    if args['schain']: suggest_chains(PDB,args)                                                                         # printing close chains for further link analysis to screen
    if type(args['debug']) is list: print_summary(PDB,output='screen')                                                  # print summary to screen
