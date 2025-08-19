#include "mcp_tools.h"
#include <iostream>
#include <jsoncpp/json/json.h>

int main(int argc, char* argv[]) {
    try {
        // Initialize MCP tools
        MCP::initializeTools();
        
        std::cout << "ExprPredictor MCP Tools Demo" << std::endl;
        std::cout << "============================" << std::endl;
        
        // Get tool registry
        MCP::ToolRegistry& registry = MCP::ToolRegistry::getInstance();
        
        // List available tools
        std::cout << "\nAvailable MCP Tools:" << std::endl;
        auto toolNames = registry.getToolNames();
        for (const auto& name : toolNames) {
            MCP::Tool* tool = registry.getTool(name);
            std::cout << "- " << name << ": " << tool->getDescription() << std::endl;
        }
        
        // Create an MCP server
        MCP::MCPServer server;
        
        // Test tools/list request
        std::cout << "\n=== Testing tools/list request ===" << std::endl;
        Json::Value listRequest;
        listRequest["method"] = "tools/list";
        
        Json::Value listResponse = server.handleRequest(listRequest);
        
        Json::StreamWriterBuilder builder;
        std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
        writer->write(listResponse, &std::cout);
        std::cout << std::endl;
        
        // Test getting schema for a specific tool
        std::cout << "\n=== Testing tool schema request ===" << std::endl;
        Json::Value schemaRequest;
        schemaRequest["method"] = "tools/get_schema";
        schemaRequest["params"]["name"] = "expr_par_get_free_pars";
        
        Json::Value schemaResponse = server.handleRequest(schemaRequest);
        writer->write(schemaResponse, &std::cout);
        std::cout << std::endl;
        
        // Test executing a tool (mock test with minimal data)
        std::cout << "\n=== Testing tool execution ===" << std::endl;
        Json::Value executeRequest;
        executeRequest["method"] = "tools/call";
        executeRequest["params"]["name"] = "expr_predictor_obj_func";
        
        // Create minimal parameters for testing
        Json::Value mockParameters;
        Json::Value mockExprPar;
        mockExprPar["maxBindingWts"] = Json::Value(Json::arrayValue);
        mockExprPar["maxBindingWts"].append(1.0);
        mockExprPar["txpEffects"] = Json::Value(Json::arrayValue);
        mockExprPar["txpEffects"].append(1.0);
        mockExprPar["repEffects"] = Json::Value(Json::arrayValue);
        mockExprPar["repEffects"].append(0.0);
        mockExprPar["basalTxp"] = 1.0;
        
        // Create minimal matrix
        Json::Value mockMatrix;
        mockMatrix["rows"] = 1;
        mockMatrix["cols"] = 1;
        mockMatrix["data"] = Json::Value(Json::arrayValue);
        Json::Value row(Json::arrayValue);
        row.append(1.0);
        mockMatrix["data"].append(row);
        mockExprPar["factorIntMat"] = mockMatrix;
        
        mockParameters["parameters"] = mockExprPar;
        executeRequest["params"]["arguments"] = mockParameters;
        
        Json::Value executeResponse = server.handleRequest(executeRequest);
        writer->write(executeResponse, &std::cout);
        std::cout << std::endl;
        
        std::cout << "\n=== MCP Tools Demo Complete ===" << std::endl;
        std::cout << "Tools are ready for use with GitHub MCP server integration." << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}