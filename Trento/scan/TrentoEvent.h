#ifndef EVENT_H
#define EVENT_H

#ifdef __CLING__
#pragma cling add_include_path("/usr/include/hdf5/serial")
#pragma cling load("/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so")
#pragma cling load("/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so")
#endif

#include <H5Cpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace trento {

class Event
{
public:
  Event();
  ~Event();

  void setEv(int i) { mEv = i; }
  void setB(double d) { mB = d; }
  void setNpart(double d) { mNpart = d; }
  void setNcoll(double d) { mNcoll = d; }
  void setmult(double d) { mMult = d; }
  void setE(std::vector<double> v) { mE = v; }

  int ev() const { return mEv; };
  double b() const { return mB; };
  double npart() const { return mNpart; };
  double ncoll() const { return mNcoll; };
  double mult() const { return mMult; };
  std::vector<double> e() const { return mE; };

  void print() const
  {
    std::cout << mEv << "\t" << mB << "\t" << mNpart << "\t" << mNcoll << "\t" << mMult;
    for (double val : mE) {
      std::cout << "\t" << val;
    }
    std::cout << std::endl;
  }

private:
  int mEv;
  double mB;
  double mNpart;
  double mNcoll;
  double mMult;
  std::vector<double> mE;
};

Event::Event()
{
}

Event::~Event()
{
}

namespace parse {

inline std::vector<Event> text(const char* infile = "data.dat", const bool printEvent = false)
{
  std::vector<Event> events;
  std::ifstream in(infile);
  if (!in.is_open()) {
    std::cerr << "ERROR: could not open input file: " << infile << std::endl;
    return events;
  }

  std::string line;
  long nRead = 0, nSkipped = 0;
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

    std::istringstream iss(line);
    int ev;
    double b, npart, ncoll, mult, e4, e5, e6, e7;
    if (iss >> ev >> b >> npart >> ncoll >> mult >> e4 >> e5 >> e6 >> e7) {
      Event evt;
      evt.setEv(ev);
      evt.setB(b);
      evt.setNpart(npart);
      evt.setNcoll(ncoll);
      evt.setmult(mult);
      evt.setE({e4, e5, e6, e7});

      events.push_back(evt);
      nRead++;
      } else {
        std::cerr << "WARNING: could not parse line, skipping: " << line << std::endl;
        nSkipped++;
    }
  }
  in.close();

  if (printEvent) {
    for (const auto& evt : events) {
      evt.print();
    }
  }

  return events;
}


inline std::vector<Event> hdf(const char* infile, const bool printEvent = false)
{
  std::vector<Event> events;
  try {
    H5::H5File file(infile, H5F_ACC_RDONLY);
    // Number of objects directly under the root group
    hsize_t nObjects = file.getNumObjs();
    events.reserve(nObjects);

    for (hsize_t i = 0; i < nObjects; ++i) {
      std::string name = file.getObjnameByIdx(i);
      H5G_obj_t type = file.getObjTypeByIdx(i);

      // We expect each event to be a dataset
      if (type != H5G_DATASET) {
          continue;
      }

      H5::DataSet dataset = file.openDataSet(name);
      Event evt;

      // Event number comes from the dataset name: "event_123"
      int ev = std::stoi(name.substr(name.find('_') + 1));
      evt.setEv(ev);

      // Helper lambda for reading a double attribute
      auto readDoubleAttr = [&](const char *attrName) {
        H5::Attribute attr = dataset.openAttribute(attrName);
        double value;
        attr.read(H5::PredType::NATIVE_DOUBLE, &value);
        return value;
      };

      // Helper lambda for reading an integer attribute
      auto readIntAttr = [&](const char *attrName) {
        H5::Attribute attr = dataset.openAttribute(attrName);
        int value;
        attr.read(H5::PredType::NATIVE_INT, &value);
        return value;
      };

      evt.setB(readDoubleAttr("b"));
      evt.setNpart(readIntAttr("npart"));
      evt.setNcoll(readIntAttr("ncoll"));
      evt.setmult(readDoubleAttr("mult"));
      evt.setE({ readDoubleAttr("e2"), readDoubleAttr("e3"), readDoubleAttr("e4"), readDoubleAttr("e5")});
      events.push_back(evt);
    }
  } catch (H5::FileIException &e) {
    std::cerr << "ERROR: could not open HDF5 file: " << infile << std::endl;
    e.printErrorStack();
    return events;
  } catch (H5::DataSetIException &e) {
    std::cerr << "ERROR: could not open an event dataset in file: " << infile << std::endl;
    e.printErrorStack();
    return events;
  } catch (H5::AttributeIException &e) {
    std::cerr << "ERROR: could not read an event attribute in file: " << infile << std::endl;
    e.printErrorStack();
    return events;
  } catch (std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return events;
  }

  if (printEvent) {
    for (const auto& evt : events) {
      evt.print();
    }
  }

  return events;
}

} // namespace parse
} // namespace trento

#endif // EVENT_H