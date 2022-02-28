#ifndef _PROPALIGNMENT_H
#define _PROPALIGNMENT_H

#include <fun4all/SubsysReco.h>

class PropAlignment: public SubsysReco
{
public:
  PropAlignment();
  virtual ~PropAlignment();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

private:
  int GetNodes(PHCompositeNode* topNode);
};

#endif
