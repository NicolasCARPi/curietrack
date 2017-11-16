/*******************************************************************************************************************
    TimeNodes Class
    Array containing all the Node s of a particular instant
    Institut Curie
    UMR - 144
    by Victor Racine
    12 july 2002
*******************************************************************************************************************/

#ifndef __time_nodes__
    #define __time_nodes__

    #include "general.h"
    #include "nodes.h"

    class TimeNodes{
    protected:
        int time;
        virtual Node* allocNode(void)=0;
    public:
        Node**	tabNodes;
        TimeNodes(int _time, int _nbNodes);
        ~TimeNodes();
        Node *operator()(int i);
        int nbNodes;
        int linkNeighborsPresent();
        int linkNeighborsBeforeAndAfter(TimeNodes & futureTimeNodes);
        int linkNeighborsBlinkBeforeAndAfter(TimeNodes & futureTimeNodes);
        int sortByIntensityAllNeighbors();
    };


#endif //__time_nodes__
