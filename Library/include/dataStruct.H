#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H


    #define STR_LENGTH_PROPINFO 30
    // some general data structure
    /**
     * @brief Information of a thermodynamic property
     * 
     */
    struct propInfo
    {
        // string shortName, longName, unit;
        // int dimension[7]; //OpenFOAM dimension
        char shortName[STR_LENGTH_PROPINFO], longName[STR_LENGTH_PROPINFO], unit[STR_LENGTH_PROPINFO];
    };

#endif