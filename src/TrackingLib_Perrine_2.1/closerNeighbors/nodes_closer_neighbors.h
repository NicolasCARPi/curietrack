/*******************************************************************************************************************
    Nodes closer neighbors Class
    It is derived from the class Nodes
    This class find the closer neighbors for each structure.
    Institut Curie
    UMR - 144
    by Victor Racine
    2004 06 21
*******************************************************************************************************************/

#ifndef __node_closer_neighbors__
#define __node_closer_neighbors__

#include "../kernel/nodes.h"

class NodeCloserNeighbors:public Node{
private:

public:
	NodeCloserNeighbors(int _time);

	void makeConnection(void);
};

#endif //__node_closer_neighbors__
