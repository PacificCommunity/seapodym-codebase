#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/XMLDocumentHandler.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

using std::cerr;
using std::cout;
using std::endl;

XERCES_CPP_NAMESPACE_USE

// class MyDOMHandler : public XMLDocumentHandler {
class MyDOMHandler : public XercesDOMParser {
   public:
    MyDOMHandler(){};
    virtual ~MyDOMHandler(){};
};

int main(int argc, char* args[]) {
    try {
        XMLPlatformUtils::Initialize();
    } catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Error during initialization! :\n" << message << "\n";
        XMLString::release(&message);
        return 1;
    }

    // XercesDOMParser* parser = new XercesDOMParser();
    MyDOMHandler* parser = new MyDOMHandler();
    parser->setValidationScheme(XercesDOMParser::Val_Always);  // optional.
    parser->setDoNamespaces(true);                             // optional

    ErrorHandler* errHandler = (ErrorHandler*)new HandlerBase();
    parser->setErrorHandler(errHandler);

    char* xmlFile = "x1.xml";

    try {
        parser->parse(xmlFile);
    } catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Exception message is: \n" << message << "\n";
        XMLString::release(&message);
        return -1;
    } catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cout << "Exception message is: \n" << message << "\n";
        XMLString::release(&message);
        return -1;
    } catch (...) {
        cout << "Unexpected Exception \n";
        return -1;
    }

    DOMDocument* doc = parser->getDocument();
    DOMElement* root = doc->getDocumentElement();
    DOMNodeIterator* iterator =
        doc->createNodeIterator(root, DOMNodeFilter::SHOW_ALL, 0, true);
    for (DOMNode* n = iterator->nextNode(); n != 0; n = iterator->nextNode()) {
        short type = n->getNodeType();
        if (type == 1) {
            char* element = XMLString::transcode(n->getNodeName());
            cout << "ELEMENT: " << strlen(element) << ':' << element << ' ';
            XMLString::release(&element);
            element = 0;

            if (n->hasAttributes()) {
                DOMNamedNodeMap* attributes = n->getAttributes();
                cout << '[';
                int size = attributes->getLength();
                for (int i = 0; i < size; i++) {
                    DOMNode* a = attributes->item(i);

                    char* attribute = XMLString::transcode(a->getNodeName());
                    cout << attribute << '=';
                    XMLString::release(&attribute);
                    attribute = 0;

                    char* value = XMLString::transcode(a->getNodeValue());
                    cout << value << ' ';
                    XMLString::release(&value);
                    value = 0;

                    a = 0;
                }

                cout << ']';

                attributes = 0;
            }

            cout << endl;
        } else if (type == 3) {
            const XMLCh* xmlch = n->getNodeValue();
            DOMText* domtext = doc->createTextNode(xmlch);
            cout << domtext->getLength() << endl;

            char* text = XMLString::transcode(domtext->getData());
            cout << "TEXT: " << strlen(text) << ':' << text << endl;
            XMLString::release(&text);
            text = 0;
            /*
                            char* text =
               XMLString::transcode(n->getNodeValue()); cout << text[0] <<
               "TEXT: " << strlen(text) << ':' << text << endl;
                            XMLString::release(&text);
                            text = 0;
            */
        }
    }

    delete parser;
    delete errHandler;

    return 0;
}
