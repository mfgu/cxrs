try:
    from numpy import *
    numx = 'numpy'
except:
    numx = 'none'

def readcol(fn, k, format=[], start=0, stop=-1, step=1,
            comment='#%', filter='', delim=None):
    """readcol reads columns from a text file

    fn is the file name

    k is a list of integers corresponding to the columns to be read.
    starting from 0

    format is a list of strings specifies the type of each column. possible values
    are 'F', 'I', 'C', 'A', corresponding to float, int, complex, and string
    respectively. If no format is specified, all columns are assumed to be float.
    if format is a scalar string, all columns are assumed to have the same format
    as that string.

    start is an integer specifying the starting line number. lines above that are
    ignored. lines counts from 0

    stop is an integer specifying the stopping line number. lines below that are
    ignored. Both start and stop may be negative, in which case, they count from
    the bottom of the file.

    step read every step lines.
    
    comment is a string specifying comment characters. if the first non-black
    character is in comment, then that line is ignored.

    filter is a string specifying a command to filter the rows. k-th field of
    the record can be referenced as c[k] in this string. e.g.,
    filter='c[0] > 10' only selects lines that have first column > 10

    delim is a string specifying the delimiter of the fields.

    The function returns a list containing the read columns. if the column is
    'F', 'I', or 'C', then the Numeric array is returned. if the column is 'A',
    then the ordinary string list is returned for that column.
    
    """
    f = open(fn, 'r')
    lines = f.readlines()
    n = len(lines)
    if (stop < 0):
        stop = n + stop
    if (start < 0):
        start = n + start
    m = len(k)
    mk = max(k)+1
    kr = range(m)
    nf = len(format)
    if (type(format) != type([]) and type(format) != type((0,))):
        format = [format]
    if (nf == 0):
        format = ['F']*m
    elif (nf != m):        
        for j in range(nf, m):
            format.append(format[-1])
    conv = []
    d = []
    for t in format:
        if t == 'F' or t == 'f':
            conv.append(float)
            if numx == 'numpy':
                d.append(zeros(n, float))
            else:
                d.append([0.0]*n)
        elif t == 'I' or t == 'i':
            conv.append(int)
            if numx == 'numpy':
                d.append(zeros(n, int))
            else:
                d.append([0]*n)
        elif t == 'C' or t == 'c':
            conv.append(complex)
            if numx == 'numpy':
                d.append(zeros(n, complex))
            else:
                d.append([0j]*n)
        elif t == 'A' or t == 'a':
            conv.append(str)
            d.append(['']*n)
        else:
            raise Exception('Illegal Format')

    p = 0
    for i in range(start, stop+1, step):
        c = lines[i].strip()
        if (len(c) == 0):
            continue
        if (c[0] in comment):
            continue
        c = c.split(delim)
        if (len(c) < mk):
            continue
        skipline = 0
        for j in kr:
            try:
                c[k[j]] = conv[j](c[k[j]])
            except:
                skipline = 1
                break
        if skipline:
            continue
        if (filter):
            if (not eval(filter)):
                continue
        for j in kr:
            d[j][p] = c[k[j]]
        p = p + 1

    f.close()

    for j in kr:
        d[j] = d[j][0:p]
        
    return d

