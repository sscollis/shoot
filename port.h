#ifdef CRAY
#define CGEEV cgeev
#else
#ifdef XLF
#define CGEEV zgeev_
#else
#define CGEEV zgeev
#endif
#endif
