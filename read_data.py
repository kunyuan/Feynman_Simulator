import numpy as np
import subprocess
import os
import sys
import StringIO

class data_list(list):
    def data(self):
        return list.__getitem__(self, 0)

    def xylabel(self):
        return list.__getitem__(self, 1)

def __index__(file, comment='#', separator='^[[:blank:]]*###*'):
    index_list = []
    name_list = []
    if os.path.exists(file):
        result=subprocess.check_output(['grep', '-b', separator, file])
        i=0
        for item in result[:-1].split('\n'):
            i+=1
            index_list.append(int(item.split(':')[0]))
            name=item.split(':')[1].strip(' \t\n\r#')
            if name=='':
                name='Block: '+str(i)+''
            name_list.append(name)
        print "Block Number: {0}".format(len(index_list))
        print "Name List: {0}".format(name_list)
        length=os.path.getsize(file)
        for i in range(0, len(index_list)-1):
            index_list[i]=(index_list[i], index_list[i+1])
        index_list[-1]=(index_list[-1], length)

    else:
        print "Data file '{0}' does not exist!".format(file)
        sys.exit()

    return index_list, name_list, length

def __parse_scope__(index_list, name_list, length,scope=':'):
    scope_bound_str=scope.split(':')
    if len(scope_bound_str)>2:
        print "Only 2 parameters are allowed in range: {}".format(scope)
        sys.exit()
    else:
        try:
            if scope_bound_str[0]!='':
                header=int(scope_bound_str[0])
            else:
                header=0

            if scope_bound_str[1]!='':
                scope_bound=index_list[header:int(scope_bound_str[1])]
                scope_name_list=name_list[header:int(scope_bound_str[1])]
            else:
                scope_bound=index_list[header:]
                scope_name_list=name_list[header:]

        except:
            print "Please check your range: {}".format(scope)
            sys.exit()

    return scope_bound,scope_name_list

def __parse_dim__(dimstr):
    dim=[]
    dim_name=[]
    strbuff=dimstr.strip(' \t\n\r')[1:].split(',')
    for i in range(0,len(strbuff)):
        elem=strbuff[i].strip(' \t\n\r').split(':')
        if len(elem)==1:
            dim.append(int(elem[0]))
            if len(strbuff)<=3:
                if i==0:
                    dim_name.append("X")
                elif i==1:
                    dim_name.append("Y")
                else:
                    dim_name.append("Z")
            else:
                dim_name.append('X'+str(i-1))

        elif len(elem)==2:
            dim.append(int(elem[1]))
            if elem[0].strip(' \t\n\r')=='':
                dim_name.append('X'+str(i-1))
            else:
                dim_name.append(elem[0])
        else:
            print "Illegal dimensional information: {}".format(dimstr)

    return dim, dim_name

def __read_one_block__(f,bound,comment='#'):
    f.seek(bound[0])
    buff=StringIO.StringIO(f.read(bound[1]-bound[0]))
    buff.readline()
    buffstr=buff.readline()
    dim, dim_name=__parse_dim__(buffstr)
    a=np.loadtxt(buff,comments=comment)
    if a.shape[1]>2:
        print "At most two numbers in each row!"
        sys.exit()
    elif a.shape[1]==2:
        a.dtype=complex
    return a.reshape(dim,order='F'),dim_name

def read_array(file, name=[], scope=':',
               comment='#', separator='^[[:blank:]]*###*'):
    '''
    Return: a dict of data blocks as {block_name: blocklist([block,dim_info]),...}
    ------------------------
        -block_name: a string contains the name of the block

        -block: an array contains data
        !!!!Attention: It always returns the last block in the file which is named after block_name
        -dim_info: a list conatins the name of each dimension, e.g. ['X','Y']

        -blocklist: a class inherits from list, including two additional methods:
         * data(): return "block"
         * xylabel(): return "dim_info"

    Parameters:
    -----------
    file : string
        The data file path and name

    name : list of string
        It contains a block name list to extract from the file

    scope : string with format 'int:int', default=':'
        The blocks range needed to read, as the index in list
    !!!! scope argument will be active only if you didn't pass any name list to name argument!!!

    comment : string, default='#'
        The comment characters in the data file

    separator : string, default='^[[:blank:]]*###*'
        The string to separate different blocks in the data file

    Usage:
    -------------
    It is recommanded to read data by passing a name list;
    if the "name" list is empty, the function will try to use "scope" argument;
    and if both "name" and "scope" are not specified, the function will return a dictionary contains the last block in the file 

    Examples:
    -------------

    block, dim_info=read_data.read_array("PATH/TO/FILE", ["Gamma","Sigma"])["Gamma"]

    is same as:

    data_block=read_data.read_array("PATH/TO/FILE", ["Gamma","Sigma"])
    block, dim_info=data_block["Gamma"]

    Or you can also use:

    block, dim_info=read_data.read_array("PATH/TO/FILE", scope="1:")["Gamma"]

    datablock=read_data.read_array("PATH/TO/FILE", scope="1:")["Gamma"]
    block = datablock.data()
    dim_info = datablock.dim_info()
    
    '''
    index_list, name_list, length = __index__(file, comment, separator)
    name_dict=dict(zip(name_list,index_list))
    #print scope_bound
    data_block={}
    f=open(file,'r')
    #print len(scope_bound)
    if len(name)==0:
        scope_bound, scope_name_list = __parse_scope__(index_list,
                                                     name_list, length, scope)
        for i in range(0,len(scope_bound)):
            block, dim_name = __read_one_block(f,scope_bound[i])
            data_block[scope_name_list[i]] = data_list([block, dim_name])
    else:
        for item in name:
            try:
                bound = name_dict[item]
            except:
                print "[Warning]The name: "+item+" doesn't exist in data file!"
                continue
            block, dim_name = __read_one_block__(f,bound)
            data_block[item] = data_list([block, dim_name])

    f.close()
    print "returned block keys list: {0}".format(data_block.keys())
    return data_block

if __name__=="__main__":
    read_array('0_0.50_1_Gam_matrix.dat',['G'])
