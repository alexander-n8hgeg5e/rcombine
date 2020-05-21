#!/usr/bin/env python3
import random
from re import sub,search
from sys import argv,exit,stderr
from math import floor,log10,ceil,floor,inf
from multiprocessing import Pool,cpu_count
from time import sleep,time
from pylib.du import ptb
from pylib.du import dd
from pprint import pprint

used_list=[]

def init_random():
    pass
    with open("/dev/urandom",'rb') as f:
        if args.verbose:
            print("# {} adding randomness ".format(time()))
        data=f.read(64)
    random.seed(data)

def get_base_values(norm='EIA E12'):
    """
    10% Standard Values (EIA E12)
    Decade multiples are available from 10 Ω through 1 MΩ
    """
    supported_norms=['EIA E12']
    if norm.lower() == supported_norms[0].lower():
        return [ 10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 68, 82 ]
    else:
        raise Exception("ERROR: Norm \""+norm+"\" not supported.\nsupported norms: "+str(supported_norms))

def gen_values():
    exp=[-2,-1,0,1,2,3,4,5]
    values=[]
    for e in exp:
        for i in get_base_values():
            values.append(i*10**e)
    return values

def print_values(values,start='[',sep=', ',end=']\n',format=False):
    lv=len(values)
    rs=start+end
    for i in range(lv):
        v=values[i]
        s0="{:.3f}".format(v)
        s=s0
        s=sub("([.][0-9]*[1-9]+)[0]*$",'\\1',s)
        s=sub("[.][0]*$",'',s)
        s=sub("000000$",'M',s)
        s=sub("000$",'K',s)
        s=sub('([1-9])([1-9])[0][0]$','\\1.\\2K',s)
        s=sub('([1-9])([1-9])[0][0]K$','\\1.\\2M',s)
        rs=rs+s+(sep if i != (lv-1) else end)
    print(rs,end='')

def res_p(resistor_values):
    """
    calc parallel resistance
    """
    r=0
    for v in resistor_values:
        r+=1/v
    return 1/r

def res_p_from_indexes(indexes):
    global values
    try:
        rp= res_p([values[i] for i in indexes])
    except IndexError as e:
        print('index error: len values={} indexes={}'.format(len(values),indexes),file=stderr)
        raise e
    return rp

def get_nearest_index(val):
    global values
    r=0
    for i in range(len(values)):
        if values[i] > val:
            above=i
            break
    if i>0 and values[i] - val > values[i-1] - val:
        return i-1
    return i

def gen_dolist(length,step,random_disable=False,fix=[]):
    do=[]
    for i in range(length):
        do.append(random.choice([-step,0,+step]))
    if len(fix) > 0:
        for i in range(length):
            for j in fix:
                if i == j:
                    do[i]=0
    if random_disable:
        keep=random.choice([set(range(0,length))-set(fix)])
        for i in range(length):
            if not i == keep:
                do[i]=0
    return do

def calc_max_step(indexes):
    lv=len(values)
    mstep=lv
    for i in indexes:
        mstep=min(mstep,lv-i-1)
        mstep=min(mstep,i)
    return mstep

def has_dupes(l):
    ll=len(l)
    for i in range(ll):
        for j in range(ll):
            if i != j:
                if l[i] == l[j]:
                    return True
    return False

def has_reserved_ones(l):
    global used_list
    if any(thing in used_list for thing in l):
        return True
    return False

def mark_used(v):
    global used_list
    used_list.append(v)

def calc_better_v2(val,indexes,delta,random_jump_lock,nodupes,exclude,worker_index):
    counter=0
    random_jump=1
    _random_jump=0
    block_jump=20
    value_index_set=set(range(len(values)))-exclude
    done=False
    while not done:
        li=len(indexes)
        for i in range(li):
            for j in value_index_set:
                if i == 0:
                    if li > 1:
                        _indexes = [j] + indexes[1:]
                    else:
                        _indexes = [j]
                elif i == li-1:
                    if li > 1:
                        _indexes = indexes[:-1] + [j]
                    else:
                        _indexes = [j]
                elif li > 2:
                    _indexes = indexes[:i] + [j] + indexes[i+1:]
                else:
                    raise Exception("ERROR: unexpected case")
                if block_jump > 0:
                    _random_jump=0
                else:
                    _random_jump = random_jump
                better = calc_better_v3(val,_indexes,delta,_random_jump, nodupes, worker_index,block_jump,fix=[i])
                if _random_jump > 0:
                    block_jump = 100
                else:
                    block_jump-=1
                if better is None:
                    pass
                else:
                    return better
        counter+=1
        random_jump=min(random_jump+1,int(len(values)/2))
        if counter%10000==0:
            init_random()

def is_better(val,new_indexes,delta):
    return abs(calc_delta(val,new_indexes)) < abs(delta)

def calc_better_v3(val,indexes,delta,random_jump,nodupes,worker_index,block_jump,fix=[]):
    lv=len(values)
    for k in range(1,random_jump+2):
        li=len(indexes)
        choices=[[ -1, +1 ]]*li
        for i in range(li):
            lci=len(choices[i])
            for j in range(lci):
                choices[i][j]*=random.randint(1,k)
        for i in range(li):
            lci=len(choices[i])
            if lci > 0 and not i in fix:
                choice=random.randint(0,lci-1)
                indexes[i] += choices[i][choice]
                if not any([i >= lv for i in  indexes]):
                    if not any([i<0 for i in indexes]):
                        if is_better(val,indexes,delta):
                            if not nodupes or not has_dupes(indexes):
                                # got new indexes at this point
                                return indexes
                indexes[i] -= choices[i][choice]
                choices[i].pop(choice)

def calc_delta(val,indexes):
    return res_p_from_indexes(indexes) - val

def calc_factor_1(base,i):
    global values
    f=values[i]/base
    ff=f-1
    fff=abs(ff)
    return fff

def calc_factor(base,indexes):
    f=0
    for i in indexes:
        f=max(calc_factor_1(base,i),f)
    return f

def dividelist(l,parts):
    ll=len(l)
    parts=min(ll,parts)
    lll=[]
    dil=[0]
    for i in range(1,parts):
        di=round(i*ll/parts)
        lll.append(l[dil[-1]:di])
        dil.append(di)
    lll.append(l[dil[-1]:])
    return lll

def print_e(e,end=" ",format=False):
    plog=abs(ceil(log10(args.precision)))+1
    formatstr='e = {:.' + str(plog) + 'f}'
    s=formatstr.format(e)
    if format:
        return s+end
    else:
        print(s,end=end)

def combine(val,max_p=3,precision=0.01,add_resistor_factor=2,random_jump_lock=200,nodupes=False):
    global values
    bases = []
    dividers = []
    for i in range(max_p-1):
        bases.append(val*(i+2))
        dividers.append(i+2)
    print("# start values: "+str(bases))
    next_base_factor = add_resistor_factor
    lb=len(bases)
    pool=Pool(cpu_count())
    for bi in range(lb):
        #print("bi = {}".format(bi))
        divider=dividers[bi]
        nearest=get_nearest_index(bases[bi])
        indexes=[nearest]*divider
        indexes_parts = dividelist(list(range(len(values))),pool._processes)
        indexes_full  = set(list(range(len(values))))
        results=[]
        r=None
        while calc_factor(bases[bi],indexes) < next_base_factor or bi==lb-1:
            delta=calc_delta(val,indexes)
            e=abs(((val+delta)/val) - 1)
            if  e <= precision:
                dupes_ok=True
                if nodupes:
                    if has_dupes(indexes) or has_reserved_ones(indexes):
                        dupes_ok=False
                if dupes_ok:
                    if nodupes:
                        for i in indexes:
                            mark_used(i)
                    #print_values([values[i] for i in indexes],sep=',  ',end=']\n')
                    return [values[i] for i in indexes],e
                else:
                    # in case precision is ok but dupes not ok
                    # need to set delta high
                    # to prevent attempts improving the the delta
                    delta=inf
            any_ready=any([i.ready() for i in results])
            if not any_ready:
                for i in range(pool._processes):
                    exclude = indexes_full - set(indexes_parts[i])
                    poolargs = [val,indexes,delta,random_jump_lock,nodupes,exclude,i]
                    r = pool.apply_async(calc_better_v2,poolargs)
                    results.append(r)
                sleep(0.01)
            done=False
            while not done:
                for res_i in range(len(results)):
                    if results[res_i].ready():
                        indexes=results[res_i].get()
                        results.pop(res_i)
                        done=True
                        break
                sleep(0.01)

def parse_args():
    from argparse import ArgumentParser
    ap=ArgumentParser()
    ap.add_argument("-p","--precision",type=float,default=0.01)
    ap.add_argument("-n","--max-resistor-count",type=int,default=3)
    ap.add_argument("-a","--add-resistor-factor",type=float,default=1,help='defines how much the resistors can differ befor another resistors is added')
    ap.add_argument("-D","--do-not-add-default-resistor-values",action="store_true",default=False)
    ap.add_argument("-R","--no-reuse-value",action="store_true",default=False, help="implicitly enabled if also -D selected")
    ap.add_argument("-v","--verbose",action="store_true",default=True)
    ap.add_argument("-o","--python-output",action="store_true",default=True,help="produce output suitable for the \"-r\" option")
    ap.add_argument("--debug",action="store_true",default=True)
    ap.add_argument(dest="target_values",nargs="+",type=str)
    ap.add_argument("--add",dest="add_values",nargs='*',type=float,default=[])
    ap.add_argument("-f","--value-file",help="file with values in a python list. The syntax is python syntax.",default=None)
    ap.add_argument("--print-values-and-exit",action="store_true")
    ap.add_argument("-u","--in-use-value-file",type=str,help="file that is pyhton \"eval-able\" and contains a list with already consumed/used values")
    global args
    args=ap.parse_args()
    if args.do_not_add_default_resistor_values:
        args.no_reuse_value = True

def resistor_string_to_float(s,enable_mult=False):
    if type(s) is str:
        if enable_mult:
            sl=s.split(' ')
            for i in range(len(sl)):
                m=search('\\d+[xX][^ ]+',sl[i])
                pos=m.start()
                d=sub('(\\d+)[xX][^ ]+','\\1',sl[i][pos:])
                d=int(d)
                obj=sub('\\d+[xX]([^ ]+).*$','\\1',sl[i][pos:])
                obj=sub('^(.*)[kK]','(\\1)*1000',obj)
                obj=sub('^(.*)[mM]','(\\1)*1000000',obj)
                i0=sl[i][:pos]
                i1=sub('\\d+[xX][^ ]+','',sl[i][pos:])
                sl[i]=i0+" ".join([obj]*d)+i1
            s="["+", ".join(sl)+"]"
        s=sub('^(.*)[kK]','(\\1)*1000',s)
        s=sub('^(.*)[mM]','(\\1)*1000000',s)
        s=eval(s)
    return float(s)

def main(): 
    global values
    parse_args()
    init_random()
    if not args.value_file is None:
        with open(args.value_file) as f:
            file_values = eval(f.read())
            file_values = [ resistor_string_to_float(v) for v in file_values ]
    else:
        file_values=[]

    if not args.do_not_add_default_resistor_values:
        add_values=[resistor_string_to_float(v) for v in args.add_values]
        values=gen_values()+args.add_values+file_values
    else:
        values= args.add_values + file_values
    values.sort()

    if not args.in_use_value_file is None:
        with open(args.in_use_value_file) as f:
            used_values = eval(f.read())
            used_values = [ resistor_string_to_float(v) for v in used_values ]
            for v in used_values:
                i=values.index(v)
                used_list.append(i)

    if args.print_values_and_exit:
        print("resistor values:")
        print("=============")
        print(values)
        print("In-use values:")
        print("=============")
        pprint([ values[i] for i in used_list ])
        exit(0)


    if args.python_output:
        print("[")

    for target_value in args.target_values:
        target_value = resistor_string_to_float(target_value,enable_mult=True)
        print("# target_value = "+str(target_value))
        vals,e=combine(target_value,precision=args.precision,max_p=args.max_resistor_count,add_resistor_factor=args.add_resistor_factor,nodupes=args.no_reuse_value)
        if args.python_output:
            print(", ".join(str(v) for v in vals)+",")
        print("######{:10.10}##{:30.30}######".format("#"*10,"#"*30))
        print("###   {:10.10}  {:30.30}   ###".format(print_e(e,format=True),str(vals)))
        print("###   {:42.42}   ###".format("final combined value = {:.3f} [Ohm]".format(res_p(vals))))
        print("######{:10.10}##{:30.30}######".format("#"*10,"#"*30))
    if args.python_output:
        print("]")
    
if __name__ == '__main__':
    main()

# vim: set foldnestmax=1 foldlevel=0 foldmethod=indent :
