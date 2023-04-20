/**
 * $Id$
 * Copyright 2005 University of Hawaii, Pelagic Fisheries Research Program
 */
#ifndef __XMLDocument_h__
#define __XMLDocument_h__

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

/// @cond TellsDoxygenToSkipThese
using std::atoi;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
// using std::atof;
using std::runtime_error;
using std::strtod;

/// @endcond

#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/XMLDocumentHandler.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/sax/SAXException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/RuntimeException.hpp>
#include <xercesc/util/XMLString.hpp>

XERCES_CPP_NAMESPACE_USE

class XMLDocument {
   public:
    XMLDocument() { initialize(); }

    XMLDocument(const XMLDocument& p) {
        cerr << "Error: " << __FILE__ << ':' << __LINE__ << '\n';
    }

    virtual ~XMLDocument() {
        if (parser != 0) {
            delete parser;
            parser = 0;
        }
        if (errHandler != 0) {
            delete errHandler;
            errHandler = 0;
        }
    }

   public:
    int read(const char* parfile) throw(runtime_error);

    string get(const string& element, const string& attribute) const;

    double getDouble(const string& element, const string& attribute) const
        throw(runtime_error) {
        double d = 0;
        string s = get(element, attribute);
        if (!s.empty()) {
            // d = atof(s.c_str());
            d = strtod(s.c_str());
        } else {
            throw runtime_error(
                "Error: XMLDocument::getDouble element=\"" + element +
                "\" does not exist\n");
        }
        return d;
    }

    int getInteger(const string& element, const string& attribute) const
        throw(runtime_error) {
        int integer = 0;
        string s = get(element, attribute);
        if (!s.empty()) {
            integer = atoi(s.c_str());
        } else {
            throw runtime_error(
                "Error: XMLDocument::getInteger  element=\"" + element +
                "\" does not exist\n");
        }
        return integer;
    }

    int write(const char* parfile) throw(runtime_error);

   private:
    /**
     * setups up xerces-c enviroment
     */
    int initialize() throw(runtime_error);
    string getValue(const DOMNode* node, const string& attribute) const;

   private:
    DOMDocument* doc;
    XercesDOMParser* parser;
    ErrorHandler* errHandler;
};
#endif
