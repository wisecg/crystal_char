CXXFLAGS += $(shell root-config --cflags) -I$(ORDIR)/Util -I$(ORDIR)/Decoders -I$(ORDIR)/IO -I$(ORDIR)/Processors -I$(ORDIR)/Management
LIBS += $(shell root-config --libs) -L$(ORDIR)/lib -lORUtil -lORDecoders -lORIO -lORProcessors -lORManagement

all: getSpectrum

getSpectrum: getSpectrum.cc
	g++ $(CXXFLAGS) -o getSpectrum getSpectrum.cc $(LIBS)

clean:
	rm getSpectrum
