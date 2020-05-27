#!/usr/bin/env python3
import random
from re import sub,search,match
from sys import argv,exit,stderr
from math import floor,log10,ceil,floor,inf
from multiprocessing import Pool,cpu_count
from time import sleep,time
from pylib.du import ptb
from pylib.du import dd
from pprint import pprint
from operator import methodcaller
from types import FunctionType
from copy import deepcopy

used_list=[]

class SingleIndex(int):
    def make_double(self,i):
        return DoubleIndex(self,i)

class DoubleIndex(tuple):

    def __new__(cls, *z, **zz):
        return tuple.__new__(cls,z,**zz)

    def __getnewargs__(self):
        return self

    def __lt__(self,i):
        lt=True
        for s in self:
            if not s < i:
                lt=False
        return lt

    def __ge__(self,i):
        ge=True
        for s in self:
            if not s >= i:
                ge=False
        return ge

class Values(list):
    """
    allow getting values by double index
    """
    def __getitem__(self,*z,**zz):
        if type(z[0])is DoubleIndex:
            l=[]
            for i in z[0]:
                if type(i) is tuple:
                    print("thing: \"{}\" type: {}".format(z[0],type(z[0])))
                    raise Exception("ERROR: wrong type, got tuple want int")
                l.append(super().__getitem__(i,**zz))
            return l
        else:
            return super().__getitem__(*z,**zz)

def test_0():
    """
    Test Values and DoubleIndex
    """
    values=gen_values()
    v=Values(values)
    il=[DoubleIndex(4,5),6]
    for i in il:
        print(v[i])
    il=[4,5,6]
    for i in il:
        print(v[i])

def test_1():
    """
    test single index \"make_double\" function
    """
    i=SingleIndex(5)
    d=i.make_double(3)
    values=gen_values()
    v=Values(values)
    print("single index:",i)
    print("double index:",d)
    dd(v[i])
    dd(v[d])
    dd(v[d[0]])
    dd(v[d[1]])

def testing():
    #import types
    #print(dir(types))
    test_1()
    pass

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

def res_s(resistor_values):
    """
    calc serial resistance
    """
    return sum_it(resistor_values)

def sum_it(thing):
    if is_iterable(thing):
        return(sum(thing))
    return thing

def res_sp(resistor_values):
    _resistor_values=[]
    for v in resistor_values:
        _resistor_values.append(res_s(v))
    return res_p(_resistor_values)

def is_iterable(thing):
    return hasattr(thing,'__iter__')

def res_sp_from_indexes(indexes):
    global values
    try:
        if type(values) is list:
            raise Exception('ERROR: wrong type')
        try:
            rp= res_sp([values[i] for i in indexes])
        except IndexError:
            print(indexes)
            raise
    except IndexError as e:
        print('INDEX ERROR: len values={} indexes={}'.format(len(values),indexes),file=stderr)
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

def calc_max_step(indexes):
    lv=len(values)
    mstep=lv
    for i in indexes:
        mstep=min(mstep,lv-i-1)
        mstep=min(mstep,i)
    return mstep

def flatten_thing(thing):
    cls=thing.__class__
    flat=cls()
    for v in thing:
        if is_iterable(v):
            flat+=cls(v)
        else:
            flat.append(v)
    return flat

def has_dupes(l):
    flat=flatten_thing(l)
    ll=len(flat)
    for i in range(ll):
        for j in range(ll):
            if i != j:
                if flat[i] == flat[j]:
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

is_better_counter=0
def is_better(val,new_indexes,delta):
    global is_better_counter
    is_better_counter+=1
    new=abs(calc_delta(val,new_indexes))
    old=abs(delta)
    if args.verbose and is_better_counter % 100000 == 0:
        print("WI = {}  CI= {:20.20} DN = {:3.2f}  DO = {:3.2f}".format(Worker_index,str(new_indexes),new,old))
    return new < old 

def indexes_not_to_high(indexes,length):
    for i in indexes:
        if type(i) is DoubleIndex:
            for j in i:
                if type(j) is list:
                    print("indexes: {}".format(indexes))
                    raise Exception("wrong type")
                if j >= length:
                    return False
        if type(i) is list:
            print("indexes: {}".format(indexes))
            raise Exception("wrong type")
        if i >= length:
            return False
    return True

def combine(val,max_p=3,precision=0.01,add_resistor_factor=2,random_jump_lock=200,nodupes=False,disable_pool=False,serial=True):
    global values
    bases = []
    dividers = []
    for i in range(max_p-1):
        bases.append(val*(i+2))
        dividers.append(i+2)
    print("# start values: "+str(bases))
    next_base_factor = add_resistor_factor
    lb=len(bases)
    _cpu_count=cpu_count()
    pool=Pool(_cpu_count)
    #pool=Pool(1)
    for bi in range(lb):
        #print("bi = {}".format(bi))
        divider=dividers[bi]
        nearest=get_nearest_index(bases[bi])
        indexes=[SingleIndex(nearest)]*divider
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
                    ##############
                    ##  return  ##
                    ##############
                    #for i in indexes:
                    #    print("index: {} val: {}".format(i,values[i]))
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
                    poolargs = [val,indexes,delta,random_jump_lock,nodupes,exclude,i,serial]
                    if not disable_pool:
                        r = pool.apply_async(calc_better_v2,poolargs)
                        results.append(r)
                sleep(0.01)
            done=False
            if disable_pool:
                done=True
                indexes=calc_better_v2(val,indexes,delta,random_jump_lock,nodupes,set([]),-1)
            while not done:
                for res_i in range(len(results)):
                    if results[res_i].ready():
                        indexes=results[res_i].get()
                        results.pop(res_i)
                        done=True
                        pool.terminate()
                        pool.join()
                        pool.close()
                        pool=Pool(_cpu_count)
                        break
                sleep(0.01)

def calc_better_v2(val,indexes,delta,random_jump_lock,nodupes,exclude,worker_index,serial):
    counter=0
    random_jump=1
    _random_jump=0
    block_jump=20
    value_index_set=set(range(len(values)))-exclude
    done=False
    global Worker_index
    Worker_index=worker_index
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
                better = calc_better_v3(val,_indexes,delta,_random_jump, nodupes, worker_index,block_jump,serial,fix=[i])
                if _random_jump > 0:
                    block_jump = 100
                else:
                    block_jump-=1
                if better is None:
                    pass
                else:
                    ##############
                    ##  return  ##
                    ##############
                    return better 

        counter+=1
        random_jump=min(random_jump+1,int(len(values)/2))
        if counter%10000==0:
            init_random()

def calc_better_v3(val,indexes,delta,random_jump,nodupes,worker_index,block_jump,serial,fix=[]):
    """
    """
    def gen_modify_func(start_index):
        def modify_func(thing,distance):
            if type(thing) is int:
                thing=SingleIndex(thing)
            # -distance , because need lower value
            if type(thing) is SingleIndex:
                ret=thing.make_double( max(start_index - distance,0) )
            elif type(thing) is DoubleIndex:
                i=random.randint(0,len(thing)-1)
                _thing=[thing[0],thing[1]]
                _thing[i]+=distance
                thing=DoubleIndex(*_thing)
                ret=thing
            return ret
        return modify_func
    lv=len(values)
    for k in range(1,random_jump+2):

        #
        #  make a list of possible modification types
        #  but do not set the extent of the modification yet
        #
        li=len(indexes)
        choices1=[ -1, +1 ]                                             # adding/subtracting
        choices=[]
        for i in range(li):
            a = choices1 
            if serial:
                a += [ gen_modify_func ]
            choices.append( a )    # 2 serial connected resistors out of one resistor

        #
        #  make plan of how much to alter the index "i"
        #
        for i in range(li):
            lci=len(choices[i])
            for j in range(lci):
                ri = random.randint(1,k)
                if serial and type(choices[i][j]) is FunctionType:
                    if choices[i][j].__name__=="gen_modify_func":
                        try:
                            cij=choices[i][j]
                            # generates modify func with start index ri ("random jump")
                            choices[i][j]=choices[i][j](ri)
                        except TypeError as e:
                            cij=choices[i][j]
                            raise e
                else:
                    choices[i][j]*=ri
        #
        # modify the indexes
        #
        for i in range(li):
            lci=len(choices[i])
            if lci > 0 and not i in fix:
                befor = indexes[i]
                choice=random.randint(0,lci-1)
                if type(choices[i][choice]) is int and type(indexes[i]) is int:
                    # add int
                    indexes[i] += choices[i][choice]
                elif type(choices[i][choice]) is FunctionType:
                    # replace with return value
                    indexes[i] = choices[i][choice](indexes[i],choice)

                #
                # check whether the modification
                # meets the requirements
                #
                if indexes_not_to_high(indexes,lv):
                    if not any([j<0 for j in indexes]):
                        if is_better(val,indexes,delta):
                            if not nodupes or not has_dupes(indexes):
                                # got new indexes at this point
                                return indexes
                indexes[i] = befor
                choices[i].pop(choice)

def check_indexes(indexes):
    raise Exception("ERRROR: check unwanted")
    fail=False
    if not indexes_not_to_high(indexes,len(values)):
        fail_reason="one of indexes to high"
        fail=True
    for i in indexes:
        if type(i) is DoubleIndex:
            # test copy
            ii=deepcopy(i)
            if not i==ii or i is ii:
                fail_reason="copy test failed {} != {}".format(i,ii)
                print(DoubleIndex(i))
                print(DoubleIndex(*i))
                fail=True
            for j in i:
                if type(j) is list:
                    fail_reason="DoubleIndex contains list"
                    fail=True
                    break
            if fail:
                break
    if fail:
        raise Exception("ERROR:index check failed. indexes: {}\nFail reason: {}".format(indexes,fail_reason))

def calc_delta(val,indexes):
    return res_sp_from_indexes(indexes) - val

def calc_factor_1(base,i):
    global values
    vi=values[i]
    f=vi/base
    ff=f-1
    fff=abs(ff)
    return fff

def calc_factor(base,indexes):
    f=0
    for i in indexes:
        if type(i) is DoubleIndex:
            l=[calc_factor_1(base,j) for j in i]
        else:
            l=[calc_factor_1(base,i)]
        f=max(l)
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

def parse_args():
    from argparse import ArgumentParser
    ap=ArgumentParser()
    ap.add_argument("-p","--precision",type=float,default=0.01,help="Factor of how much the combined value is allowed to diff from the target value.")
    ap.add_argument("-n","--max-resistor-count",type=int,default=3,help="How many branches/resistors in \"parallel\" you want to use at max to combine to one target value.")
    ap.add_argument("-a","--add-resistor-factor",type=float,default=1,help='Factor of how much the values in the set of the combination of resistors can diff from each other befor the length of the set is enlarged and another place for a resistor is added to the search.')
    ap.add_argument("-D","--do-not-add-default-resistor-values",action="store_true",default=False,help="Implicitly selects \"--no-reuse-value\" and no default values are added.")
    ap.add_argument("-R","--no-reuse-value",action="store_true",default=False, help="Implicitly enabled if also -D selected. No values will be used more than once. This is for use with a custom value list in contrast to use a default list.")
    ap.add_argument("-v","--verbose",action="store_true",default=True)
    ap.add_argument("--testing",action="store_true",default=False)
    ap.add_argument("-o","--python-output",action="store_true",default=True,help="produce output suitable for the \"-r\" option")
    ap.add_argument("--debug",action="store_true",default=True)
    ap.add_argument("-P","--disable_pool",action="store_true",default=False)
    ap.add_argument("-s","--serial",action="store_true",default=False,help="additionally use serial resistor combinations")
    ap.add_argument(dest="target_values",nargs="+",type=str,help='possible syntax := \"5k 10K 2x1k 2000.5 50.3k 5x2.3M 2.3m\". Whereby \"2x\" means two times the following thing (no space). EXAMPLE CALL:\"./rcombine.py -DR -f ./my_available_resistors -n2 -p0.01 -u ./found_resistors_combinations_that_will_not_be_used_in_the_search 2x5.0k 5x10.0K\"')
    ap.add_argument("--add",dest="add_values",nargs='*',type=float,default=[])
    ap.add_argument("-f","--value-file",help="These are the values that you want to use to be combined to the target values/value . Specify a file with values in a python list. The syntax is python syntax.",default=None)
    ap.add_argument("--print-values-and-exit",action="store_true",help="Print the values that are available to be combined to the target values. The output depends on the the other options that affect the available values. Use this if you want to check the list of values to find out if the the program understands your input as desired.")
    ap.add_argument("-u","--in-use-value-file",type=str,help="File that is python \"eval-able\" and contains a list with already consumed/used values. Paste the output of this program with the -o option into this file. The values from this file will not be used in the search of combinations.")
    global args
    args=ap.parse_args()
    if args.do_not_add_default_resistor_values:
        args.no_reuse_value = True

def resistor_string_sub_symbols(s):
    s=sub('([ ]|^)([^ ]+)[kK]','\\1(\\2)*1000',s)
    s=sub('([ ]|^)([^ ]+)[mM]','\\1(\\2)*1000000',s)
    return s

def resistor_string_to_float(s,enable_mult=False):
    if type(s) is str:
        if enable_mult:
            parts=s.split(' ')
            sl=[]
            for part in parts:
                if match('^\\d+[xX][^ ]+$',part):
                    d=sub('^(\\d+)[xX][^ ]+$','\\1',part)
                    d=int(d)
                    obj=sub('^\\d+[xX]([^ ]+)$','\\1',part)
                    sl+=([obj]*d)
                else:
                    sl.append(part)
            for i in range(len(sl)):
                sl[i]=float(eval(resistor_string_sub_symbols(sl[i])))
            s=sl
        else:
            s=float(eval(resistor_string_sub_symbols(s)))
    else:
        s=float(s)
    return s

def prepare_target_values(values,enable_mult=True):
    _values=[]
    for v in values:
        _values.append( resistor_string_to_float(  v, enable_mult=True ))
    values=[]
    for v in _values:
        if hasattr(v,'__iter__'):
            values+=prepare_target_values(v)
        else:
            values.append(v)
    return values

def main(): 
    global values
    parse_args()
    if args.testing:
        testing()
        exit()
    init_random()
    if not args.value_file is None:
        with open(args.value_file) as f:
            _file_values = eval(f.read())
            file_values=[]
            for v in _file_values:
                float_or_floats=resistor_string_to_float(v)
                if hasattr(float_or_floats,'__iter__'):
                    if not type(float_or_floats) is list:
                        float_or_floats=list(float_or_floats)
                        file_values+=float_or_floats
                else:
                    file_values.append(float_or_floats)
    else:
        file_values=[]

    if not args.do_not_add_default_resistor_values:
        add_values=[resistor_string_to_float(v) for v in args.add_values]
        values=gen_values()+args.add_values+file_values
    else:
        values = args.add_values + file_values
    values.sort()
    values=Values(values)

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

    target_values=prepare_target_values(args.target_values,enable_mult=True)
    for target_value in target_values:
        print("# target_value = "+str(target_value))
        vals,e=combine  (
                        target_value,
                        precision=args.precision,
                        max_p=args.max_resistor_count,
                        add_resistor_factor=args.add_resistor_factor,
                        nodupes=args.no_reuse_value,
                        disable_pool=args.disable_pool,
                        serial=args.serial,
                        )
        if args.python_output:
            flat=flatten_thing(vals)
            print(", ".join(str(v) for v in flat)+",")
        print_result(vals,e)
    if args.python_output:
        print("]")

def print_result(vals,e):
    lines=[]
    lines.append( [print_e(e,format=True),str(vals)])
    lines.append( [("final combined value = {:.3f} [Ohm]".format(res_sp(vals)))])
    width=0
    for line in lines:
        line_width=0
        for thing in line:
            line_width+=len(thing)
        line_width+=len(line)-1
        width=max(line_width,width)
    print("#"*(width+12))
    for line in lines:
        fs='###   {:'+str(width)+"."+str(width)+'}   ###'
        print(fs.format(" ".join(line)))
    print("#"*(width+12))

    
if __name__ == '__main__':
    main()

# vim: set foldnestmax=1 foldlevel=0 foldmethod=indent :
