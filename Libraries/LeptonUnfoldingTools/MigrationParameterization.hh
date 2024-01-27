#ifndef MigrationParameterization_hh
#define MigrationParameterization_hh

class TCanvas;
class TF1;

class MigrationParameterization {

public:
  MigrationParameterization() { }
  virtual ~MigrationParameterization() { }

  virtual TCanvas* DrawParameterCanvas() { return nullptr; }

  virtual TF1* MakeFunction(double E) const = 0;
  virtual bool IsPdf() const = 0;

  virtual double UnderflowProbability(double) const { return 0.0; }
  virtual double OverflowProbability(double) const { return 0.0; }

};

#endif
