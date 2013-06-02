// ***************************************************************************
// bamtools_resolve.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 23 June 2011
// ---------------------------------------------------------------------------
// Resolves paired-end reads (marking the IsProperPair flag as needed).
// ***************************************************************************

#ifndef BAMTOOLS_ILLUMINA_H
#define BAMTOOLS_ILLUMINA_H

#include "bamtools_tool.h"

namespace BamTools {

class IlluminaTool : public AbstractTool {

    public:
        IlluminaTool(void);
        ~IlluminaTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]);

    private:
        struct IlluminaSettings;
        IlluminaSettings* m_settings;

        struct IlluminaToolPrivate;
        IlluminaToolPrivate* m_impl;

        struct ReadNamesFileReader;
        struct ReadNamesFileWriter;
        struct StatsFileReader;
        struct StatsFileWriter;
};

} // namespace BamTools

#endif // BAMTOOLS_ILLUMINA_H
