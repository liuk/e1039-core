#ifndef PHG4ParameterContainerInterface__H
#define PHG4ParameterContainerInterface__H

#include <map>
#include <string>

class PHCompositeNode;
class PHG4Parameters;
class PHG4ParametersContainer;

//! This class is deprecated, please use PHG4Parameter instead.
//! See https://github.com/sPHENIX-Collaboration/coresoftware/pull/405
class PHG4ParameterContainerInterface
{
 public:
  PHG4ParameterContainerInterface(const std::string &name);
  virtual ~PHG4ParameterContainerInterface();

  void set_name(const std::string &name);
  virtual void  SetDefaultParameters() = 0;

 // Get/Set parameters from macro
  void set_double_param(const int id, const std::string &name, const double dval);
  double get_double_param(const int id, const std::string &name) const;
  void set_int_param(const int id, const std::string &name, const int ival);
  int get_int_param(const int id, const std::string &name) const;
  void set_string_param(const int id, const std::string &name, const std::string &sval);
  std::string get_string_param(const int id, const std::string &name) const;

  void UpdateParametersWithMacro();
  void CreateInitialize(const int detid);
  void SaveToNodeTree(PHCompositeNode *runNode, const std::string &nodename);
  void PutOnParNode(PHCompositeNode *parNode, const std::string &nodename);
  int ExistDetid(const int detid) const;

 protected:
  void set_default_double_param( const std::string &name, const double dval);
  void set_default_int_param( const std::string &name, const int ival);
  void set_default_string_param( const std::string &name, const std::string &sval);
  void InitializeParameters();
  const PHG4ParametersContainer *GetParamsContainer() {return paramscontainer;}
  PHG4ParametersContainer *GetParamsContainerModify() {return paramscontainer;}
  const PHG4Parameters *GetDefaultParameters() {return defaultparams;}

 private:
  PHG4ParametersContainer *paramscontainer;
  PHG4Parameters *defaultparams;
  std::map<int, PHG4Parameters *> macroparams;
};

#endif
