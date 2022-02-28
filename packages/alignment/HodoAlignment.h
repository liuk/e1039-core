#ifndef _HODOALIGNMENT_H
#define _HODOALIGNMENT_H

#include <fun4all/SubsysReco.h>

class HodoAlignment: public SubsysReco
{
public:
  HodoAlignment();
  virtual ~HodoAlignment();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

private:
  int GetNodes(PHCompositeNode* topNode);
};

#endif
