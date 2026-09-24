import sys, types
class _Dummy:
    def __init__(self,*a,**k): pass
    def __call__(self,*a,**k): return _Dummy()
    def __getattr__(self,n): return _Dummy()
def _mod(name):
    m=types.ModuleType(name)
    def ga(n):
        if n.startswith('__'): raise AttributeError(n)
        return type(n,(_Dummy,),{})
    m.__getattr__=ga
    return m
for n in ['tkinter','tkinter.ttk','tkinter.messagebox','tkinter.filedialog','tkinter.font','tkinter.scrolledtext','tkinter.simpledialog','matplotlib.backends.backend_tkagg','_tkinter']:
    sys.modules[n]=_mod(n)
t=sys.modules['tkinter']; t.ttk=sys.modules['tkinter.ttk']; t.messagebox=sys.modules['tkinter.messagebox']; t.filedialog=sys.modules['tkinter.filedialog']
for c in ['END','BOTH','X','Y','LEFT','RIGHT','TOP','BOTTOM','W','E','N','S','NW','NE','SW','SE','CENTER','DISABLED','NORMAL','RAISED','SUNKEN','FLAT','GROOVE','RIDGE','HORIZONTAL','VERTICAL','WORD','NONE','INSERT','ACTIVE','YES','NO','SINGLE','EXTENDED','BROWSE','MULTIPLE']:
    setattr(t,c,c.lower())
