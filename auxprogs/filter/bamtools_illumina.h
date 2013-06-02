// ***************************************************************************
// bamtools_illumina 
// Tonatiuh Pena-Centeno
// ---------------------------------------------------------------------------
// Last modified:  2-June-2013
// ---------------------------------------------------------------------------
// Filters out paired-end reads (so far Illumina technology is considered).
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

