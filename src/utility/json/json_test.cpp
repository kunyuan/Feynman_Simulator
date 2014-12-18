//
//  json_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/16/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "json.h"
#include <iostream>
#include <fstream>
using namespace std;

void TestJson()
{

    Json::Value root; // will contains the root value after parsing.
    Json::Reader reader;
    std::ifstream test("../test.json", ios::in);
    test >> root;
    //    bool parsingSuccessful = reader.parse(test, root);
    //    if (!parsingSuccessful) {
    //        // report to the user the failure and their locations in the document.
    //        std::cout << "Failed to parse configuration\n"
    //                  << reader.getFormattedErrorMessages();
    //        return;
    //    }

    // Get the value of the member of root named 'encoding', return 'UTF-8' if there is no
    // such member.
    Json::Value interaction = root["Interaction"];
    cout << interaction[1].asFloat() << endl;
    // Get the value of the member of root named 'encoding', return a 'null' value if
    // there is no such member.
    //    const Json::Value plugins = root["plug-ins"];
    //    for (int index = 0; index < plugins.size(); ++index) // Iterates over the sequence elements.
    //        loadPlugIn(plugins[index].asString());
    //
    //    setIndentLength(root["indent"].get("length", 3).asInt());
    //    setIndentUseSpace(root["indent"].get("use_space", true).asBool());
    //
    //    // ...
    //    // At application shutdown to make the new configuration document:
    //    // Since Json::Value has implicit constructor for all value types, it is not
    //    // necessary to explicitly construct the Json::Value object:
    //    root["encoding"] = getCurrentEncoding();
    //    root["indent"]["length"] = getCurrentIndentLength();
    //    root["indent"]["use_space"] = getCurrentIndentUseSpace();

    Json::StyledWriter writer;
    // Make a new JSON document for the configuration. Preserve original comments.
    std::string outputConfig = writer.write(root);

    // You can also use streams.  This will put the contents of any JSON
    // stream at a particular sub-value, if you'd like.
    //    std::cin >> root["subtree"];
}