/**
 * $Id$
 * Copyright@2005 University of Hawaii, Pelagic Fisheries Research Program
 */
#include "XMLDocument.h"

int XMLDocument::initialize() throw (runtime_error) {
	try {
		XMLPlatformUtils::Initialize();
	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		//cout << "Error during initialization! :\n" << message << "\n";
		throw runtime_error(message);
		XMLString::release(&message);
		return -1; ///< mycomment
	}

	parser = new XercesDOMParser();
	parser->setValidationScheme(XercesDOMParser::Val_Always);
	parser->setDoNamespaces(true);

	errHandler = (ErrorHandler*) new HandlerBase();
	parser->setErrorHandler(errHandler);

	return 1;
}

int XMLDocument::read(const char* xmlFile) throw (runtime_error) {
	try {
		parser->parse(xmlFile);
	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		//cout << "Exception message is: \n" << message << "\n";
		throw runtime_error(message);
		XMLString::release(&message);
		return -1;
	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		//cout << "Exception message is: \n" << message << "\n";
		throw runtime_error(message);
		XMLString::release(&message);
		return -1;
	} catch (const SAXException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		//cout << "Exception message is: \n" << message << "\n";
		throw runtime_error(message);
		XMLString::release(&message);
		return -1;
	}

	return 1;
}

string XMLDocument::get(const string& element, const string& attribute) const {
	string value;
	
	DOMDocument *doc = parser->getDocument();
	DOMElement *root = doc->getDocumentElement();

	DOMNodeIterator *iterator = doc->createNodeIterator(root, DOMNodeFilter::SHOW_ALL, 0, true);
	for (DOMNode *node = iterator->nextNode(); value.empty() && node != 0; node = iterator->nextNode()) {
		char* nodeelement = XMLString::transcode(node->getNodeName());

		if (strcmp(nodeelement, element.c_str()) == 0) { 
			value = getValue(node, attribute);
		}

		XMLString::release(&nodeelement);
		nodeelement = 0;
	}

	return value;
}

string XMLDocument::getValue(const DOMNode* node, const string& attribute) const {
	string value;

	const short type = node->getNodeType();
	if (type == 1) {
		DOMNamedNodeMap *attributes = node->getAttributes();
		if (attributes != 0) {
			const int size = attributes->getLength();
			for (int i = 0; value.empty() && i < size; i++) {
				DOMNode* a = attributes->item(i);
				char* nodeattribute = XMLString::transcode(a->getNodeName());

				if (strcmp(nodeattribute, attribute.c_str()) == 0) { 
					char* v = XMLString::transcode(a->getNodeValue());
					value = v;
					XMLString::release(&v);
					v = 0;
				}

				XMLString::release(&nodeattribute);
				nodeattribute = 0;
				a = 0;
			}
		}
		attributes = 0;
	}

	return value;
}

int write(const char* parfile) throw (runtime_error) {
	return 0;
}
