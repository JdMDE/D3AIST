#ifndef DEB_H
#define DEB_H

// WARNING: 
// When using multiline macros, remember that the line-sepparation character \ MUST be the last character of its line. 
#define ERRPROG(x) \
{\
    std::cerr << "Error from ";\
    if (progname!="")\
        std::cerr << progname << ": ";\
    std::cerr << x\
    exit(1);\
}

#define WARNPROG(x) \
{\
    std::cerr << "Warning from ";\
    if (progname!="")\
        std::cerr << progname << ": ";\
    std::cerr << x\
    std::cerr.flush();\
}
#endif
