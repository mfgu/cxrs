from sys import version_info
if version_info.major == 2:
    from Tkinter import *
    from ScrolledText import *
    from FileDialog import *
    from popen2 import *
else:
    from tkinter import *
    from tkinter.scrolledtext import *
    from tkinter.filedialog import *
    from subprocess import *
import os, glob, select, string
from textio import readcol
import re

# this class removes the unix hidden dir from the file dialog
class MyFileDialog(FileDialog):    
    def filter_command(self, event=None):
        dir, pat = self.get_filter()
        try:
            names = os.listdir(dir)
        except os.error:
            self.master.bell()
            return
        self.directory = dir
        self.set_filter(dir, pat)
        names.sort()
        subdirs = [os.pardir]
        matchingfiles = []
        for name in names:
            fullname = os.path.join(dir, name)
            if os.path.isdir(fullname):
                subdirs.append(name)
            elif fnmatch.fnmatch(name, pat):
                matchingfiles.append(name)
        self.dirs.delete(0, END)
        for name in subdirs:
            if (name[0] != '.' or name == '.' or name == '..'):
                self.dirs.insert(END, name)
        self.files.delete(0, END)
        for name in matchingfiles:
            self.files.insert(END, name)
        head, tail = os.path.split(self.get_selection())
        if tail == os.curdir: tail = ''
        self.set_selection(tail)

# a scrolled listbox
class ScrolledListbox(Listbox):
    def __init__(self, master, **kw):
        self.frame = Frame(master)
        self.vbar = Scrollbar(self.frame, orient=VERTICAL)
        self.vbar.pack(side=RIGHT, fill=Y)
        Listbox.__init__(self, self.frame, **kw)
        self.pack(side=LEFT, fill=BOTH, expand=1)
        self['yscrollcommand'] = self.vbar.set
        self.vbar['command'] = self.yview
        
        # Copy geometry methods of self.frame -- hack!
        methods = list(Pack.__dict__.keys())
        methods = methods + list(Grid.__dict__.keys())
        methods = methods + list(Place.__dict__.keys())

        for m in methods:
            if m[0] != '_' and m != 'config' and m != 'configure':
                setattr(self, m, getattr(self.frame, m))
                
class PixList:
    def __init__(self, master, top):
        self.topwin = top
        self.default = [0,1,3,4,7,11,13,14,16,17,19,
                        20,21,22,23,24,25,27,29,30,31];
        f = Frame(master)
        f.pack(side=LEFT)
        lab = Label(f, text='Pixels', width=8)
        lab.grid(row=0, column=0)
        self.ref = Entry(f, width=8)
        self.ref.grid(row=0, column=1)
        self.ref.insert(END, self.default[0])
        b0 = Listbox(f, selectmode=MULTIPLE, exportselection=0,
                     width=8, height=16)
        b0.grid(row=1, column=0)
        b1 = Listbox(f, selectmode=MULTIPLE, exportselection=0,
                     width=8, height=16)
        b1.grid(row=1, column=1)
        for i in range(16):
            b0.insert(END, 'Pix %02d'%i)
            b1.insert(END, 'Pix %02d'%(i+16))
        self.list0 = b0
        self.list1 = b1
        if version_info.major == 2:
            w = 5
        else:
            w = 6
        b = Button(f, text='All', width=w, command=self.check_list0)
        b.grid(row=2, column=0)
        b = Button(f, text='None', width=w, command=self.uncheck_list0)
        b.grid(row=3, column=0)
        b = Button(f, text='Reset', width=w, command=self.reset_list0)
        b.grid(row=4, column=0)
        b = Button(f, text='All', width=w, command=self.check_list1)
        b.grid(row=2, column=1)
        b = Button(f, text='None', width=w, command=self.uncheck_list1)
        b.grid(row=3, column=1)
        b = Button(f, text='Reset', width=w, command=self.reset_list1)
        b.grid(row=4, column=1)

    def check_list0(self):
        self.list0.select_set(0, END)
        
    def check_list1(self):
        self.list1.select_set(0, END)

    def uncheck_list0(self):
        self.list0.select_clear(0, END)

    def uncheck_list1(self):
        self.list1.select_clear(0, END)

    def reset_list0(self):
        self.uncheck_list0()
        for i in self.default:
            if i < 16:
                self.list0.select_set(i)
                
    def reset_list1(self):
        self.uncheck_list1()
        for i in self.default:
            if i >= 16:
                self.list1.select_set(i-16)

    def reset(self):
        self.reset_list0()
        self.reset_list1()
        
    def ref_pix(self):
        i = int(self.ref.get())
        p = self.pixels()
        if not i in p:
            i = p[0]
            self.ref.delete(0, END)
            self.ref.insert(END, '%d'%i)
            s = 'Illegal reference pixel, reset to %d\n'%i
            self.topwin.update_log(s)
        return '%d'%i

    def pixels(self):
        p0 = list(map(int, self.list0.curselection()))
        p1 = list(map(int, self.list1.curselection()))
        for i in range(len(p1)):
            p1[i] += 16
        return p0+p1

class DirList:
    def __init__(self, master, top):
        self.topwin = top
        f = Frame(master)
        f.pack(side=LEFT)
        f0 = Frame(f)
        f0.pack(side=LEFT)
        f1 = Frame(f)
        f1.pack(side=LEFT)
        f2 = Frame(f1)
        f2.pack()
        b = Label(f2, text='HistD', width=4)
        b.pack(side=LEFT)
        b = Entry(f2, width=11)        
        b.pack()
        self.ref = b
        d1 = ScrolledListbox(f1, selectmode=MULTIPLE, exportselection=0,
                             width=13, height=21)
        d1.pack()
        self.hlist = d1
        b = Label(f0, text='FList', width=14)
        b.pack()
        d0 = ScrolledListbox(f0, selectmode=MULTIPLE, exportselection=0,
                             width=15, height=16)
        d0.pack()
        self.dlist = d0
        f1 = Frame(f0)
        f1.pack()
        if version_info.major == 2:
            w = 5
        else:
            w = 5
        b = Button(f1, text='Update', width=w, command=self.update_dir)
        b.grid(row=0, column=0)
        b = Button(f1, text='Reset', width=w, command=self.reset)
        b.grid(row=0, column=1)
        b = Button(f1, text='AllF', width=w, command=self.all_list)
        b.grid(row=1, column=0)
        b = Button(f1, text='AllH', width=w, command=self.all_hist)
        b.grid(row=1, column=1)
        b = Button(f1, text='NoneF', width=w, command=self.none_list)
        b.grid(row=2, column=0)
        b = Button(f1, text='NoneH', width=w, command=self.none_hist)
        b.grid(row=2, column=1)

    def update_dir(self):
        hdir = self.topwin.opts.get_opt('HDir')
        fns = glob.glob('%s[0-9][0-9]'%hdir)
        cfn = self.dlist.get(0, END)
        for i in range(len(cfn)-1, -1, -1):
            if not cfn[i] in fns:
                self.dlist.delete(i)
        cfn = self.dlist.get(0, END)
        fns.sort()
        for a in fns:
            if not a in cfn:
                self.dlist.insert(END, a)
                self.dlist.select_set(END)

        cfn = self.dlist.get(0, END)
        s = list(map(int, self.dlist.curselection()))
        h = []
        ih = []
        for i in range(len(cfn)):
            a = cfn[i]
            b = readcol(a, [2], format=['A'])
            b = b[0]
            h += b
            ih += [i in s]*len(b)
            
        ch = self.hlist.get(0, END)
        for i in range(len(ch)-1, -1, -1):
            if not ch[i][:-3] in cfn:
                self.hlist.delete(i)                
        ch = self.hlist.get(0, END)
        for i in range(len(ch)):
            try:
                j = cfn.index(ch[i][:-3])
                if not j in s:
                    self.hlist.select_clear(i)
            except:
                pass
        h.sort()
        for i in range(len(h)):
            a = h[i]
            if not a in ch:
                self.hlist.insert(END, a)
                if (ih[i]):
                    self.hlist.select_set(END)
                else:
                    self.hlist.select_clear(END)
        a = self.hlist.get(0)
        self.ref.delete(0, END)
        if a:
            self.ref.insert(END, a)
        
    def reset(self):
        self.update_dir()
        self.all_list()
        self.all_hist()
        
    def all_list(self):
        self.dlist.select_set(0, END)

    def all_hist(self):
        self.hlist.select_set(0, END)

    def none_list(self):
        self.dlist.select_clear(0, END)

    def none_hist(self):
        self.hlist.select_clear(0, END)
        
    def ref_dir(self):        
        d = self.ref.get()
        return d
    
    def dirs(self):
        cf = self.dlist.get(0, END)
        ch = self.hlist.get(0, END)
        sh = list(map(int, self.hlist.curselection()))
        h = [ch[i] for i in sh]
        tmin = [0.0]*len(h)
        tmax = [0.0]*len(h)
        for a in cf:
            r = readcol(a, [0,1,2], format=['F','F','A'])
            for i in range(len(h)):
                try:
                    j = r[2].index(h[i])                    
                    tmin[i] = r[0][j]
                    tmax[i] = r[1][j]
                except:
                    pass
        r = [h,tmin,tmax]
        return r
    
class EvtList:
    def __init__(self, master, top):
        f = Frame(master)
        f.pack(side=LEFT)
        lab = Label(f, text='EvtList', width=5).pack()
        f0 = Frame(f)
        f0.pack(side=BOTTOM)
        if version_info.major == 2:
            w0 = 5
            w1 = 3
        else:
            w0 = 6
            w1 = 4
        b = Button(f0, text='Update', width=w0, command=self.update_evt)
        b.pack(side=LEFT)
        b = Button(f0, text='All', width=w1, command=self.select_all)
        b.pack(side=LEFT)
        b = Button(f0, text='None', width=w1, command=self.clear_all)
        b.pack(side=LEFT)        
        b = ScrolledListbox(f, selectmode=MULTIPLE, exportselection=0,
                            width=22, height=19)
        b.pack()
        self.list = b
        self.topwin = top
    
    def update_evt(self):
        self.topwin.set_defaults()
        self.topwin.pixlist.reset()
        fns = self.evtfiles
        cfn = self.list.get(0, END)
        for i in range(len(cfn)-1,-1,-1):
            self.list.delete(i)
        for i in range(len(fns)):
            a = fns[i]
            self.list.insert(END, a)
            self.list.select_set(END)

    def select_all(self):
        self.list.select_set(0, END)

    def clear_all(self):
        self.list.select_clear(0, END)
        
    def evt_files(self):
        s = list(map(int, self.list.curselection()))
        return s

class Opts:
    def __init__(self, master, top, side=LEFT):
        self.topwin = top
        f = Frame(master)
        f.pack(side=side)
        self.opts_name = ['HDir', 'VBin', 'EBin', 'DRange',
                          'Scale', 'Control', 'IRange', 'PRange',
                          'FPrep', 'FFit', 'CDir', 'FMerge', 'FLine', 'FPoly',
                          'NPoly', 'FPlot', 'DMerge', 'FDCorr', 'NDCorr',
                          'ADCorr', 'CEvent', 'CRef']
        self.opts_default0 = ['d', '1e-3, 0.0, 15.0, 1',
                              '0.0, 0.0, 1.5e4', '0.5, 10',
                              '0.0', '0.0', '0', '0.5, 10, 6, 5.0',
                              'Ln,Kn', '', '', '',
                              'lines.txt', 'cpc', '5', '',
                              '', '', '500,0,0', '.dc,1,0,1', '.ce,0',
                              '00t00']
        n = len(self.opts_name)
        if n%2 == 1:
            self.opts_name.append('Dummy')
            self.opts_default0.append('')
            n = n+1
        self.opts_default = self.opts_default0[:]
        self.opts_val = []        
        n2 = int(n/2)
        for i in range(n2):
            Label(f, text=self.opts_name[i]).grid(row=i, column=0)
            e = Entry(f, width=32)
            e.grid(row=i, column=1)
            self.opts_val.append(e)
        for i in range(n2, n):
            Label(f, text=self.opts_name[i]).grid(row=i-n2, column=2)
            e = Entry(f, width=32)
            e.grid(row=i-n2, column=3)
            self.opts_val.append(e)
        
    def reset(self):
        for i in range(len(self.opts_name)):
            self.opts_val[i].delete(0, END)
            self.opts_val[i].insert(END, self.opts_default[i])
            
    def get_opt(self, optn):
        try:
            i = self.opts_name.index(optn)
            v = self.opts_val[i].get()
        except:
            raise Exception('%s not found'%optn)

        v = v.strip()
        if not v:
            raise Exception('%s: No value set'%optn)
        # validate the option values.
        if optn == 'HDir':
            v = v.strip()
            if len(v) > 12:
                raise Exception('%s: HDir too long, must <= 12 char'%optn)
        if optn == 'VBin':
            try:
                x = v.split(',')
                if (len(x) != 4):
                    raise Exception('Invalid VBin')
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))
        if optn == 'EBin':
            try:
                x = v.split(',')
                if (len(x) < 3):
                    raise Exception('Invalid EBin')
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))
        elif optn == 'DRange':
            try:
                x = v.split(',')
                if (len(x) <= 1):
                    raise Exception('Invalid DRange')
                x = [float(a) for a in x]
                for j in range(len(x)-1):
                    if x[j+1] <= x[j]:
                        raise Exception('Invalid DRange')
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))            
        elif optn == 'PRange':
            try:
                x = v.split(',')
                if (len(x) != 4):
                    raise Exception('Invalid PRange')
                x = [float(a) for a in x]         
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))                   
        elif optn == 'Scale':
            try:
                x = v.split(',')
                if (len(x) < 1):
                    raise Exception('Invalid Scale')
                x = [float(a) for a in x]
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))
        elif optn == 'Control':
            try:
                x = v.split(',')
                if (len(x) < 1):
                    raise Exception('Invalid DRange')
                x = [float(a) for a in x]
                for j in range(len(x)-1):
                    if x[j+1] <= x[j]:
                        raise Exception('Invalid DRange')
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))     
        elif optn == 'NPoly':
            try:
                xv = v.split(',');
                x = int(xv[0])
                if x < 0:
                    raise Exception('Invalid NPoly')
                if (len(xv) > 1):
                    x = int(xv[1])
                    if (x < 0 or x > 1):
                        raise Exception('Invalid NPoly')
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))
        elif optn == 'ADCorr':
            try:
                x = v.split(',')
                if (len(x) != 4):
                    raise Exception('Invalid ADCorr')
                x = [int(a) for a in x[1:]]
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))
        elif optn == 'CDir':
            try:
                v = v.strip()
            except Exception as e:
                raise Exception('%s: %s'%(optn, str(e)))
        elif optn == 'CRef':
            try:
                v = v.strip()
                if len(v) > 8:
                    raise Exception('%s: CRef too long, must <= 8 char'%opn)
            except Exception as e:
                raise Exception('%s: %s'%(optn,str(e)))                    
        return v
    
class TopDir:
    def __init__(self, master, top, evtlistf):
        self.master = master
        f = Frame(master)
        f.pack()
        self.entry = ent = Entry(f, width=120)
        ent.pack(side=LEFT)
        ent.bind("<Return>", self.chdir)
        if version_info.major == 2:
            w = 4
        else:
            w = 5
        bro = Button(f, text='Load', width=w, command=self.choose_dir)
        bro.pack(side=LEFT)
        sav = Button(f, text='Save', width=w, command=self.save_file)
        sav.pack(side=LEFT)
        self.topwin = top
        self.evtlistf = evtlistf
        self.entry.insert(END, os.path.join(os.getcwd(), self.evtlistf))
        self.fullevtlistf = None
        
        Button(f, text='Cancel', width=w, command=top.cancel).pack(side=LEFT)
        Button(f, text='Quit', width=w, command=master.quit).pack(side=LEFT)
        
    def chdir(self, event=None):        
        cwd = os.getcwd()
        ofn = self.evtlistf
        fn = self.entry.get()
        if fn == self.fullevtlistf:
            return
        if fn:
            fn = os.path.expanduser(fn)
            fn = os.path.abspath(fn)
            if os.path.isfile(fn):
                (d, self.evtlistf) = os.path.split(fn)
            elif os.path.isdir(fn):
                d = fn
            else:
                d = cwd
            if not d:
                d = cwd
            if d != cwd:
                try:
                    os.chdir(d)            
                    s = os.path.join(d, self.evtlistf)
                    self.entry.delete(0, END)
                    self.entry.insert(END, s)
                    s = 'CWD: evtlist file changed to %s.\n'%s
                    self.topwin.update_log(s)
                    self.topwin.reset()
                except:
                    self.evtlistf = ofn
                    self.entry.delete(0, END)
                    self.entry.insert(END, os.path.join(os.getcwd(),self.evtlistf))
                    self.topwin.update_log('CWD: cannot change to %s.\n'%d)
            else:
                s = os.path.join(d, self.evtlistf)
                self.entry.delete(0, END)
                self.entry.insert(END, s)
                s = 'CWD: evtlist file changed to %s.\n'%s
                self.topwin.update_log(s)
                self.topwin.reset()
            self.fullevtlistf = self.entry.get()
            
    def choose_dir(self):
        d = MyFileDialog(self.master, title='Loading Event List File')
        fn = d.go(key='evtload')
        if (fn):
            self.entry.delete(0, END)
            self.entry.insert(END, fn)
            self.chdir()

    def save_file(self):
        r = MyFileDialog(self.master, title='Saving Event List File')
        ofn = r.go(key='evtsave')
        if not ofn:
            return
        fn = self.entry.get()
        try:
            f = open(fn, 'r')
            d = f.read()
            f.close()
            d = re.sub(r'\r', '\n', d)
            d = d.split('\n')
            for i in range(len(d)):
                a = d[i].strip()
                if not a:
                    continue
                if a[0] != '#':
                    break
            d = d[i:]
        except:
            d = ['# evtfile    nhists  phase0  phase1  filter\n']
        
        try:
            f = open(ofn, 'w')
            ipx = self.topwin.pixlist.pixels()
            if len(ipx):
                f.write('#Pixel %d'%ipx[0])
                for i in range(1, len(ipx)):
                    f.write(' ,%d'%ipx[i])
                f.write('\n')
            for i in range(len(self.topwin.opts.opts_name)):
                a = self.topwin.opts.opts_name[i]
                try:
                    b = self.topwin.opts.get_opt(a)
                except:
                    b = ''
                if b:
                    s = '#%s %s\n'%(a,b)
                    f.write(s)

            for a in d:
                f.write('%s\n'%a)
            f.close()
            s = 'evtlist file %s saved\n'%ofn
            self.topwin.update_log(s)
        except Exception as e:
            s = 'Error in saving %s\n'%ofn
            self.topwin.update_log(s)
        
class CXRS:
    def __init__(self, master, evtlistf):
        self.root = master
        self.topdir = TopDir(master, self, evtlistf)
        f1 = Frame(master)
        f1.pack()
        self.logw = ScrolledText(f1, width=76)
        self.logw.insert(END, 'XRS Calibration Program\n')
        self.logw.pack(side=RIGHT)      
        f0 = Frame(f1)
        f0.pack(side=LEFT)
        f2 = Frame(f0)
        f2.pack()
        f = Frame(f0)
        f.pack()
        self.opts = Opts(f, self)
        n = len(self.opts.opts_name)
        self.logw.configure(height=int(30+n*0.82))
        f = Frame(master)
        f.pack()
        if version_info.major == 2:
            w = 5
        else:
            w = 5
        Button(f, text='Reset', width=w, command=self.reset0).pack(side=LEFT)
        Button(f, text='MkHist', width=w, command=self.mkhist).pack(side=LEFT)
        Button(f, text='Calib', width=w, command=self.calib).pack(side=LEFT)
        Button(f, text='AddSp', width=w, command=self.addsp).pack(side=LEFT)
        Button(f, text='PrepL', width=w, command=self.prepln).pack(side=LEFT)
        Button(f, text='FitSp', width=w, command=self.fitsp).pack(side=LEFT)
        Button(f, text='MergeL', width=w, command=self.mergeln).pack(side=LEFT)
        Button(f, text='TabV', width=w, command=self.tabvolts).pack(side=LEFT)
        Button(f, text='CPoly', width=w, command=self.cpoly).pack(side=LEFT)
        Button(f, text='MergeP', width=w, command=self.mergepoly).pack(side=LEFT)
        Button(f, text='APoly', width=w, command=self.apoly).pack(side=LEFT)
        Button(f, text='PltSp', width=w, command=self.pltsp).pack(side=LEFT)
        Button(f, text='DCorr', width=w, command=self.dcorr).pack(side=LEFT)
        Button(f, text='MergeD', width=w, command=self.merged).pack(side=LEFT)
        Button(f, text='ApplyD', width=w, command=self.applyd).pack(side=LEFT)
        Button(f, text='CEvent', width=w, command=self.cevent).pack(side=LEFT)
        Button(f, text='MEvent', width=w, command=self.mevent).pack(side=LEFT)
        
        self.pixlist = PixList(f2, self)
        self.evtlist = EvtList(f2, self)
        self.dirlist = DirList(f2, self)

        self.topdir.chdir()
        self.cmdpid = None

    def reset(self):
        self.evtlist.update_evt()
        self.opts.reset()
        self.dirlist.reset()

    def reset0(self):
        self.topdir.fullevtlistf = None
        self.topdir.chdir()
        
    def set_defaults(self):
        fn = self.topdir.evtlistf
        try:
            f = open(fn, 'r')
            d = f.read()
            f.close()
            if not d:
                raise Exception('Empty file %s'%fn)
            d = re.sub(r'\\[ \t]*[\r\n]', ' ', d)
            d = re.sub(r'\r', '\n', d)
            d = re.sub(r'[ \t]+', ' ', d)
            f = open(fn+'.m2u', 'w')
            f.write(d)
            if (d[-1] != '\n'):
                f.write('\n')
            f.close()
            d = d.split('\n')
            pix = []
            evtf = []
            self.opts.opts_default = self.opts.opts_default0[:]
            for a in d:
                a = a.strip()
                if not a:
                    continue
                if (a[0] == '#'):
                    a += ' '*20
                    if (a[1:6].lower() == 'pixel'):
                        x = a[6:].strip().split(',')
                        for y in x:
                            try:
                                pix.append(int(y))
                            except:
                                pass
                    else:
                        for i in range(len(self.opts.opts_name)):
                            optn = self.opts.opts_name[i].lower()
                            n = len(optn)
                            if a[1:n+1].lower() == optn:
                                self.opts.opts_default[i] = a[n+1:].strip()
                else:
                    x = a.split()
                    if len(x) < 5:
                        s = 'insufficient columns in evtlist, %d'%len(x)
                        raise Exception(s)
                    if not os.path.isfile(x[0]):
                        raise Exception('cannot find evt file %s'%x[0])
                    evtf.append(x[0])
                        
            if (len(pix)):
                self.pixlist.default = pix
            self.evtlist.evtfiles = evtf
            
            i = self.opts.opts_name.index('FFit')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i] = '%s,0.01'%self.opts.opts_default[0]
            i = self.opts.opts_name.index('FMerge')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i] = '%s,100'%self.opts.opts_default[0]
            i = self.opts.opts_name.index('FPlot')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i] = '%s.csp'%self.opts.opts_default[0]
            i = self.opts.opts_name.index('FDCorr')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i] = '%s.dcc'%self.opts.opts_default[0]
            i = self.opts.opts_name.index('DMerge')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i] = '%s.dcc0'%self.opts.opts_default[0]
            i = self.opts.opts_name.index('CDir')
            if not self.opts.opts_default[i]:
                self.opts.opts_default[i]='NULL'
                
        except Exception as e:
            s = os.path.join(os.getcwd(), fn)
            s = 'cannot read evtlist %s properly\n'%s            
            self.update_log(s)
            s = 'Error: %s\n'%str(e)
            self.update_log(s)
            self.evtlist.evtfiles = []
            
    def strlist2str(self, s):
        a = '['
        for i in range(len(s)):
            if (i == 0):
                a += '"%s"'%s[i]
            else:
                a += ',"%s"'%s[i]
        a += ']'
        return a

    def cancel(self):
        try:
            if (self.cmdpid != None):
                a = os.system('kill -9 %s'%(self.cmdpid))
        except:
            pass
        
    def executecmd(self, c, ofn):
        self.topdir.chdir()
        if version_info.major == 2:
            a = Popen4(c, 0)
            cio = [a.fromchild]
        else:
            a = Popen([c], shell=True, bufsize=0, stderr=PIPE, stdin=PIPE, stdout=PIPE)
            cio = [a.stdout, a.stderr]
        self.cmdpid = a.pid
        f = open(ofn, 'w')
        try:
            p = select.poll()
            for x in cio:
                p.register(x)
        except:
            p = None

        wt = 50
        while (1):
            fc = 0
            r = a.poll()
            if (p):
                b = p.poll(wt)
                for d in b:
                    if (d[1] &
                        (select.POLLNVAL | select.POLLHUP | select.POLLERR)):
                        fc = 1
                        break
                    if (d[1] & select.POLLIN):
                        s = os.read(d[0], 1024)
                        if version_info.major == 3:
                            s = str(s, 'utf-8')
                        f.write(s)
                        f.flush()
                        self.update_log(s)
            else:
                b = select.select(cio, [], [], wt*0.001)
                b = b[0]
                for d in b:
                    s = d.readline()
                    f.write(s)
                    f.flush()
                    self.update_log(s)
            if (r != -1 and r != None):
                break
            if (fc):
                r = a.poll()
                break
            self.root.update()
        self.update_log('Done: %s, ret=%s\n'%(c, r))
        f.close()
        self.cmdpid = None

    def splitlines(self, s):
        s1 = ''
        inquote = 0
        mc = 50
        nc = 0
        for i in range(len(s)):
            s1 += s[i]
            nc += 1
            if s[i] == '"':
                inquote = (not inquote)
            if (not inquote) and (s[i] == ','):
                if nc >= mc:
                    s1 += '\n'
                    nc = 0
        return s1
    
    def mkhist(self):
        try:
            evtf = self.evtlist.evt_files()
            if (len(evtf) == 0):
                raise Exception('No event files')
            evtlistf = self.topdir.evtlistf
            hdir = self.opts.get_opt('HDir').strip()
            a = self.opts.get_opt('VBin')
            a = a.split(',')
            dv = a[0]
            vmin = a[1]
            vmax = a[2]
            md = a[3]
            a = self.opts.get_opt('FDCorr')
            a = a.split(',')
            fdc = a[0].strip()
            nd = 0
            if len(a) > 1:
                nd = int(a[1])
            a = self.opts.get_opt('CEvent')
            a = a.split(',')
            iph = a[1].strip()
            sfn = 'mkhist_%s.sl'%hdir
            ofn = 'mkhist_%s.log'%hdir            
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return

        c = 'filter_args.use_ph = %s'%iph
        f.write('%s;\n'%c)
        if nd:
            a = (evtf, evtlistf+'.m2u', hdir, dv, vmin, vmax, md, fdc)
            c = 'splitevts(%s, "%s", "%s", %s, %s, %s, NULL, %s, "%s")'%a
        else:
            a = (evtf, evtlistf+'.m2u', hdir, dv, vmin, vmax, md)
            c = 'splitevts(%s, "%s", "%s", %s, %s, %s, NULL, %s, NULL)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script mkhist_%s'%hdir
        self.executecmd(c, ofn)
        self.dirlist.update_dir()

    def mevent(self):
        try:
            evtf = self.evtlist.evt_files()
            if (len(evtf) == 0):
                raise Exception('No event files')            
            evtlistf = self.topdir.evtlistf
            hdir = self.opts.get_opt('HDir').strip()
            a = self.opts.get_opt('CEvent')
            a = a.split(',')
            try:
                iph = a[1].strip()
                oefn = a[2].strip()
            except:
                raise Exception('CEvent option must have 3rd arg for output')
            sfn = 'mevent_%s.sl'%hdir
            ofn = 'mevent_%s.log'%hdir            
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return 
        
        c = 'filter_args.use_ph = %s'%iph
        f.write('%s;\n'%c)
        a = (evtf, evtlistf+'.m2u', oefn)
        c = 'mevent(%s, "%s", "%s")'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script mevent_%s'%hdir
        self.executecmd(c, ofn)
        self.dirlist.update_dir()
        
    def calib(self):
        try:
            sr = self.pixlist.ref_pix()
            ipx = self.pixlist.pixels()
            hdir = self.opts.get_opt('HDir').strip()
            dh = self.dirlist.dirs()
            dh = dh[0]
            if (len(dh) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh)
            dr = self.dirlist.ref_dir()
            if (not dr):
                dr = 'NULL'
            else:
                dr = '"%s"'%dr
            xr = self.opts.get_opt('DRange')
            pr = self.opts.get_opt('PRange')
            xri = self.opts.get_opt('IRange')
            xrr = self.opts.get_opt('Scale')
            wsq = self.opts.get_opt('Control')
            sfn = 'calib_%s.sl'%hdir
            ofn = 'calib_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        
        a = (sr, ipx, dh, dr, xr, pr, xri, xrr, wsq)
        a = '%s, %s, %s, %s, [%s], [%s], [%s], [%s], [%s]'%a
        c = 'calib1d(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script calib_%s'%hdir
        self.executecmd(c, ofn)

    def addsp(self):
        try:
            cfn = self.opts.get_opt('HDir').strip()
            dh = self.dirlist.dirs()
            if (len(dh[0]) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh[0])
            dc = self.opts.get_opt('CDir').strip()
            dc = self.strlist2str(dc.split(','))
            ipx = self.pixlist.pixels()
            rc = self.opts.get_opt('CRef').strip()
            rc = self.strlist2str(rc.split(','))
            sfn = 'addsp_%s.sl'%cfn
            ofn = 'addsp_%s.log'%cfn
            f = self.open_script(sfn)
            rd = self.dirlist.ref_dir()
            if (not rd):
                rd = 'NULL'
            else:
                rd = '"%s"'%rd
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return

        a = (cfn, dh, ipx, dc, rc, rd)
        a = '"%s", %s, %s, %s, %s, %s'%a
        c = 'addspec(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close();
        c = 'isis-script addsp_%s'%cfn
        self.executecmd(c, ofn)

    def tabvolts(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            fline = self.opts.get_opt('FLine')
            fpoly = self.opts.get_opt('FPoly')
            fpoly = fpoly.split(',')
            fpoly = fpoly[0]
            dh = self.dirlist.dirs()
            ipx = self.pixlist.pixels()
            dh = dh[0]
            if (len(dh) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh)
            cdir = self.opts.get_opt('CDir')
            cdir = self.strlist2str(cdir.split(','))
            cref = self.opts.get_opt('CRef')
            cref = self.strlist2str(cref.split(','))
            dr = self.dirlist.ref_dir()
            if (not dr):
                dr = 'NULL'
            else:
                dr = '"%s"'%dr
            sfn = 'tabvolts_%s.sl'%hdir
            ofn = 'tabvolts_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        
        a = (fpoly, fline, dh, ipx, cdir, cref, dr)
        a = '"%s", "%s", %s, %s, %s, %s, %s'%a
        c = 'tabvolts(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script tabvolts_%s'%hdir
        self.executecmd(c, ofn)

    def cpoly(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            fpoly = self.opts.get_opt('FPoly')
            fpoly = fpoly.split(',')
            fpoly = fpoly[0]
            npv = self.opts.get_opt('NPoly').split(',')
            np = npv[0]
            if (len(npv) > 1) :
                md = npv[1]
            else:
                md = '0'
            dh = self.dirlist.dirs()
            dh = dh[0]
            if (len(dh) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh)
            sfn = 'cpoly_%s.sl'%hdir
            ofn = 'cpoly_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        
        a = (fpoly, np, dh, md)
        a = '"%s", %s, %s, %s'%a
        c = 'calibpoly(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script cpoly_%s'%hdir
        self.executecmd(c, ofn)

    def mergepoly(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            fpoly = self.opts.get_opt('FPoly')
            fpoly = fpoly.split(',')
            fpoly = self.strlist2str(fpoly)
            dh = self.dirlist.dirs()
            dh = dh[0]
            if (len(dh) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh)
            sfn = 'mergepoly_%s.sl'%hdir
            ofn = 'mergepoly_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return

        a = (fpoly, dh)
        a = '%s, %s'%a
        c = 'mergepoly(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script mergepoly_%s'%hdir
        self.executecmd(c, ofn)
        
    def prepln(self):
        try:
            hdir = self.opts.get_opt('HDir')
            lnf = self.opts.get_opt('FPrep').split(',')
            lnf = self.strlist2str(lnf)
            pr = self.opts.get_opt('PRange')
            pr = pr.split(',')
            x0 = pr[0]
            x1 = pr[1]
            sfn = 'prepln_%s.sl'%hdir
            ofn = 'prepln_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        
        a = (hdir, lnf, x0, x1);
        a = '"%s", %s, %s, %s'%a
        c = 'prepln(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script prepln_%s'%hdir
        self.executecmd(c, ofn)

    def fitsp(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            lnf = self.opts.get_opt('FPrep').split(',')
            lnf = self.strlist2str(lnf)
            a = self.opts.get_opt('FFit')
            b = a.split(',')
            if (len(b) >= 2):
                a = b[0]
                c = b[1].strip().split("&");
                try:
                    w = float(c[0])
                except:
                    raise Exception('Invalid FFit width in FitSp')
                alpha = -1
                nb = 1
                if len(c) > 1:
                    if c[1].strip() == '':
                        alpha = -1.0
                    else:
                        try:
                            alpha = float(c[1])
                        except:
                            raise Exception('Invalid FFit alpha in FitSp')
                    if len(c) > 2:
                        if c[2].strip() == '':
                            nb = 1
                        else:
                            try:
                                nb = int(c[2])
                            except:
                                raise Exception('Invalid FFit nb in FitSp')
                if (len(b) >= 3):
                    e = b[2].strip()
                else:
                    e = ''
            elif (len(b) == 1):
                w = 0.0
            else:
                raise Exception('Invalid FFit in FitSp')
            a = a.strip()
            pr = self.opts.get_opt('PRange')
            pr = pr.split(',')
            x0 = pr[0]
            x1 = pr[1]
            if (type(w) != str and w <= -1e-3):
                if (1+w == 1):
                    a = hdir
                w = self.pixlist.ref_pix()
                x1 = self.dirlist.ref_dir().strip()[-5:]
                x1 = '"%s"'%x1
                pfn = self.opts.get_opt('FPoly')+'.cpc'
                pfn = '"%s"'%pfn
            else:
                pfn = 'NULL'
            sfn = 'fitsp_%s.sl'%hdir
            ofn = 'fitsp_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        
        a = (hdir, lnf, a, e, x0, x1, w, pfn, alpha, nb)
        a = '"%s", %s, "%s", "%s", %s, %s, %s, %s, %s, %s'%a
        c = 'fitsp(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script fitsp_%s'%hdir
        self.executecmd(c, ofn)

    def mergeln(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            fline = self.opts.get_opt('FLine')
            a = self.opts.get_opt('FMerge')
            a = a.split(',')
            for i in range(len(a)):
                try:                
                    cnts = int(a[i])
                    lab = a[:i]
                    ext = a[i+1:]
                    break
                except:
                    pass
            lab = [x.strip() for x in lab]
            if (len(lab) == 0):
                raise Exception('Needs line files in FMerge Option')
            mfn = []
            for i in range(len(lab)):
                mfn.append(lab[i])
                try:
                    mfn[-1] = mfn[-1] + ext[i].strip()
                except:
                    pass
                a = mfn[-1]
                if not os.path.isfile(a+'.lno'):
                    raise Exception('Line list for %s does not exist'%a)
                
            mfn = self.strlist2str(mfn)
            lab = self.strlist2str(lab)
            sfn = 'mergeln_%s.sl'%hdir
            ofn = 'mergeln_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
            
        a = (fline, mfn, lab, cnts)
        a = '"%s", %s, %s, %s'%a
        c = 'mergeln(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script mergeln_%s'%hdir
        self.executecmd(c, ofn)

    def apoly(self):
        try:
            ipx = self.pixlist.pixels()
            hdir = self.opts.get_opt('HDir').strip()
            dh = self.dirlist.dirs()
            dh = dh[0]
            if (len(dh) == 0):
                raise Exception('No hist dir to process')
            dh = self.strlist2str(dh)
            a = self.opts.get_opt('EBin')
            a = a.split(',')
            dv = a[0]
            vmin = a[1]
            vmax = a[2]
            fpoly = self.opts.get_opt('FPoly')+'.cpc'
            sfn = 'apoly_%s.sl'%hdir
            ofn = 'apoly_%s.log'%hdir            
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
            
        a = (fpoly, dh, hdir+'.csp', dv, vmin, vmax, ipx)
        a = '"%s", %s, "%s", %s, %s, %s, %s'%a
        c = 'applygain(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script apoly_%s'%hdir
        self.executecmd(c, ofn)

    def pltsp(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            fn = self.opts.get_opt('FPlot')
            sfn = 'pltsp_%s.sl'%hdir
            ofn = 'pltsp_%s.log'%hdir            
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        a = (fn,)
        a = '"%s"'%a
        c = 'pltsp(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script pltsp_%s'%hdir
        self.executecmd(c, ofn)

    def dcorr(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            evtf = self.evtlist.evt_files()
            if (len(evtf) == 0):
                raise Exception('No event files')
            evtlistf = self.topdir.evtlistf
            nd = self.opts.get_opt('NDCorr')
            a = nd.split(',')
            if len(a) != 3:
                raise Exception('NDCorr must contain 3 ints in dcorr')
            fp = self.opts.get_opt('FPoly')
            a = self.opts.get_opt('CEvent')
            a = a.split(',')
            iph = a[1].strip()
            sfn = 'dcorr_%s.sl'%hdir
            ofn = 'dcorr_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        c = 'filter_args.use_ph = %s'%iph
        f.write('%s;\n'%c)
        a = (hdir, evtf, evtlistf+'.m2u', fp+'.cpc', nd)
        a = '"%s", %s, "%s", "%s", [%s]'%a
        c = 'dcorr(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script dcorr_%s'%hdir
        self.executecmd(c, ofn)

    def merged(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            evtf = self.evtlist.evt_files()
            mfn = self.opts.get_opt('DMerge')
            mfn = mfn.split(',')
            if (len(mfn) == 0):
                raise Exception('Needs drift files in DMerge Option')
            mfn = self.strlist2str(mfn)
            fp = self.opts.get_opt('FDCorr')
            fp = fp.split(',')[0]
            sfn = 'merged_%s.sl'%hdir
            ofn = 'merged_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        a = (evtf, fp, mfn)
        a = '%s, "%s", %s'%a
        c = 'merged(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script merged_%s'%hdir
        self.executecmd(c, ofn)
        
    def applyd(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            evtf = self.evtlist.evt_files()
            evtlistf = self.topdir.evtlistf
            fp = self.opts.get_opt('FDCorr')
            a = self.opts.get_opt('ADCorr')
            a = a.split(',')
            edc = a[0]
            idc = a[1]
            ift = a[2]
            iwt = a[3]
            a = self.opts.get_opt('CEvent')
            a = a.split(',')
            iph = a[1].strip()
            sfn = 'applyd_%s.sl'%hdir
            ofn = 'applyd_%s.log'%hdir
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return
        c = 'filter_args.use_ph = %s'%iph
        f.write('%s;\n'%c)
        a = (evtf, evtlistf, fp, idc, ift, iwt, edc)
        a = '%s, "%s", "%s", %s, %s, %s, "%s"'%a
        c = 'applyd(%s)'%a
        c = self.splitlines(c)
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script applyd_%s'%hdir
        self.executecmd(c, ofn)

    def cevent(self):
        try:
            hdir = self.opts.get_opt('HDir').strip()
            evtf = self.evtlist.evt_files()
            evtlistf = self.topdir.evtlistf
            fpoly = self.opts.get_opt('FPoly') + '.cpc'
            a = self.opts.get_opt('CEvent')
            a = a.split(',')
            ef = a[0].strip()
            iph = a[1].strip()
            sfn = 'cevent_%s.sl'%hdir
            ofn = 'cevent_%s.log'%hdir            
            f = self.open_script(sfn)
        except Exception as e:
            self.update_log('%s\n'%str(e))
            return

        c = 'filter_args.use_ph = %s'%iph
        f.write('%s;\n'%c)
        a = (evtf, evtlistf, fpoly, ef)
        a = '%s, "%s", "%s", "%s"'%a
        c = 'cevent(%s)'%a
        f.write('%s;\n'%c)
        f.close()
        c = 'isis-script cevent_%s'%hdir
        self.executecmd(c, ofn)
        
    def open_script(self, fn):
        f = open(fn, 'w')
        f.write('() = evalfile("xrscalib");\n')
        return f
    
    def update_log(self, s):
        if version_info.major == 2:
            s = s.translate(string.maketrans('\r', ' '))
        else:
            s = s.translate(str.maketrans('\r', ' '))
        self.logw.insert(END, s)
        self.logw.see(END)
        self.logw.update_idletasks()

root = Tk()
root.title('XRS Calibration')
try:
    evtf = sys.argv[1]
except:
    evtf = 'evtlist.txt'
    
x = CXRS(root, evtf)
    
root.mainloop()
