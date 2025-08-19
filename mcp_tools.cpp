#include "mcp_tools.h"
#include <iostream>
#include <memory>
#include <stdexcept>

namespace MCP {

// ToolRegistry implementation
ToolRegistry& ToolRegistry::getInstance() {
    static ToolRegistry instance;
    return instance;
}

void ToolRegistry::registerTool(std::unique_ptr<Tool> tool) {
    std::string name = tool->getName();
    tools_[name] = std::move(tool);
}

Tool* ToolRegistry::getTool(const std::string& name) {
    auto it = tools_.find(name);
    return (it != tools_.end()) ? it->second.get() : nullptr;
}

std::vector<std::string> ToolRegistry::getToolNames() const {
    std::vector<std::string> names;
    for (const auto& pair : tools_) {
        names.push_back(pair.first);
    }
    return names;
}

Json::Value ToolRegistry::getToolsSchema() const {
    Json::Value schema(Json::arrayValue);
    for (const auto& pair : tools_) {
        schema.append(pair.second->getSchema());
    }
    return schema;
}

// Serialization utilities
namespace Serialization {

Json::Value vectorToJson(const std::vector<double>& vec) {
    Json::Value result(Json::arrayValue);
    for (double val : vec) {
        result.append(val);
    }
    return result;
}

Json::Value matrixToJson(const Matrix& matrix) {
    Json::Value result;
    result["rows"] = matrix.nRows();
    result["cols"] = matrix.nCols();
    Json::Value data(Json::arrayValue);
    
    for (int i = 0; i < matrix.nRows(); i++) {
        Json::Value row(Json::arrayValue);
        for (int j = 0; j < matrix.nCols(); j++) {
            row.append(matrix.getElement(i, j));
        }
        data.append(row);
    }
    result["data"] = data;
    return result;
}

Json::Value intMatrixToJson(const IntMatrix& matrix) {
    Json::Value result;
    result["rows"] = matrix.nRows();
    result["cols"] = matrix.nCols();
    Json::Value data(Json::arrayValue);
    
    for (int i = 0; i < matrix.nRows(); i++) {
        Json::Value row(Json::arrayValue);
        for (int j = 0; j < matrix.nCols(); j++) {
            row.append(matrix(i, j));
        }
        data.append(row);
    }
    result["data"] = data;
    return result;
}

std::vector<double> jsonToVector(const Json::Value& json) {
    std::vector<double> result;
    for (const auto& val : json) {
        result.push_back(val.asDouble());
    }
    return result;
}

Matrix jsonToMatrix(const Json::Value& json) {
    int rows = json["rows"].asInt();
    int cols = json["cols"].asInt();
    Matrix result(rows, cols);
    
    const Json::Value& data = json["data"];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.setElement(i, j, data[i][j].asDouble());
        }
    }
    return result;
}

IntMatrix jsonToIntMatrix(const Json::Value& json) {
    int rows = json["rows"].asInt();
    int cols = json["cols"].asInt();
    IntMatrix result(rows, cols);
    
    const Json::Value& data = json["data"];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.setElement(i, j, data[i][j].asInt());
        }
    }
    return result;
}

Json::Value exprParToJson(const ExprPar& par) {
    Json::Value result;
    result["maxBindingWts"] = vectorToJson(par.maxBindingWts);
    result["factorIntMat"] = matrixToJson(par.factorIntMat);
    result["txpEffects"] = vectorToJson(par.txpEffects);
    result["repEffects"] = vectorToJson(par.repEffects);
    result["basalTxp"] = par.basalTxp;
    result["nFactors"] = par.nFactors();
    return result;
}

ExprPar jsonToExprPar(const Json::Value& json) {
    std::vector<double> maxBindingWts = jsonToVector(json["maxBindingWts"]);
    Matrix factorIntMat = jsonToMatrix(json["factorIntMat"]);
    std::vector<double> txpEffects = jsonToVector(json["txpEffects"]);
    std::vector<double> repEffects = jsonToVector(json["repEffects"]);
    double basalTxp = json["basalTxp"].asDouble();
    
    return ExprPar(maxBindingWts, factorIntMat, txpEffects, repEffects, basalTxp);
}

} // namespace Serialization

// ExprPredictorObjFuncTool implementation (simplest to start with)
ExprPredictorObjFuncTool::ExprPredictorObjFuncTool() 
    : Tool("expr_predictor_obj_func", "Compute objective function value for ExprPredictor with given parameters") {}

Json::Value ExprPredictorObjFuncTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    
    Json::Value inputSchema;
    inputSchema["type"] = "object";
    inputSchema["properties"]["parameters"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["parameters"]["type"] = "object";
    inputSchema["properties"]["parameters"]["description"] = "ExprPar parameters object";
    inputSchema["required"] = Json::Value(Json::arrayValue);
    inputSchema["required"].append("parameters");
    
    schema["inputSchema"] = inputSchema;
    return schema;
}

Json::Value ExprPredictorObjFuncTool::execute(const Json::Value& parameters) {
    try {
        if (!parameters.isMember("parameters")) {
            throw std::runtime_error("Missing 'parameters' field");
        }
        
        ExprPar par = Serialization::jsonToExprPar(parameters["parameters"]);
        
        // Note: This is a simplified implementation
        // In practice, you'd need an ExprPredictor instance with proper initialization
        // For now, we'll return a mock result
        Json::Value result;
        result["success"] = true;
        result["value"] = 0.0; // Placeholder value
        result["message"] = "Objective function computation completed (mock implementation)";
        
        return result;
        
    } catch (const std::exception& e) {
        Json::Value error;
        error["success"] = false;
        error["error"] = e.what();
        return error;
    }
}

// ExprParGetFreeParsTool implementation
ExprParGetFreeParsTool::ExprParGetFreeParsTool()
    : Tool("expr_par_get_free_pars", "Extract free parameters from ExprPar object") {}

Json::Value ExprParGetFreeParsTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    
    Json::Value inputSchema;
    inputSchema["type"] = "object";
    inputSchema["properties"]["parameters"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["parameters"]["type"] = "object";
    inputSchema["properties"]["parameters"]["description"] = "ExprPar parameters object";
    
    inputSchema["properties"]["coopMat"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["coopMat"]["type"] = "object";
    inputSchema["properties"]["coopMat"]["description"] = "Cooperativity matrix";
    
    inputSchema["properties"]["actIndicators"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["actIndicators"]["type"] = "array";
    inputSchema["properties"]["actIndicators"]["description"] = "Activator indicators";
    
    inputSchema["properties"]["repIndicators"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["repIndicators"]["type"] = "array";
    inputSchema["properties"]["repIndicators"]["description"] = "Repressor indicators";
    
    Json::Value required(Json::arrayValue);
    required.append("parameters");
    required.append("coopMat");
    required.append("actIndicators");
    required.append("repIndicators");
    inputSchema["required"] = required;
    
    schema["inputSchema"] = inputSchema;
    return schema;
}

Json::Value ExprParGetFreeParsTool::execute(const Json::Value& parameters) {
    try {
        if (!parameters.isMember("parameters") || !parameters.isMember("coopMat") ||
            !parameters.isMember("actIndicators") || !parameters.isMember("repIndicators")) {
            throw std::runtime_error("Missing required fields");
        }
        
        ExprPar par = Serialization::jsonToExprPar(parameters["parameters"]);
        IntMatrix coopMat = Serialization::jsonToIntMatrix(parameters["coopMat"]);
        
        // Convert boolean arrays from JSON
        std::vector<bool> actIndicators;
        for (const auto& val : parameters["actIndicators"]) {
            actIndicators.push_back(val.asBool());
        }
        
        std::vector<bool> repIndicators;
        for (const auto& val : parameters["repIndicators"]) {
            repIndicators.push_back(val.asBool());
        }
        
        std::vector<double> freePars;
        par.getFreePars(freePars, coopMat, actIndicators, repIndicators);
        
        Json::Value result;
        result["success"] = true;
        result["freePars"] = Serialization::vectorToJson(freePars);
        result["message"] = "Free parameters extracted successfully";
        
        return result;
        
    } catch (const std::exception& e) {
        Json::Value error;
        error["success"] = false;
        error["error"] = e.what();
        return error;
    }
}

// Placeholder implementations for other tools
ExprPredictorTrainTool::ExprPredictorTrainTool()
    : Tool("expr_predictor_train", "Train ExprPredictor model with given data") {}

Json::Value ExprPredictorTrainTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    // Add more detailed schema later
    return schema;
}

Json::Value ExprPredictorTrainTool::execute(const Json::Value& parameters) {
    Json::Value result;
    result["success"] = false;
    result["error"] = "Not implemented yet";
    return result;
}

ExprPredictorPredictTool::ExprPredictorPredictTool()
    : Tool("expr_predictor_predict", "Predict expression values using trained ExprPredictor") {}

Json::Value ExprPredictorPredictTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    return schema;
}

Json::Value ExprPredictorPredictTool::execute(const Json::Value& parameters) {
    Json::Value result;
    result["success"] = false;
    result["error"] = "Not implemented yet";
    return result;
}

ExprParLoadTool::ExprParLoadTool()
    : Tool("expr_par_load", "Load ExprPar parameters from file") {}

Json::Value ExprParLoadTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    
    Json::Value inputSchema;
    inputSchema["type"] = "object";
    inputSchema["properties"]["filename"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["filename"]["type"] = "string";
    inputSchema["properties"]["filename"]["description"] = "Path to parameter file";
    
    inputSchema["properties"]["coopMat"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["coopMat"]["type"] = "object";
    inputSchema["properties"]["coopMat"]["description"] = "Cooperativity matrix";
    
    inputSchema["properties"]["repIndicators"] = Json::Value(Json::objectValue);
    inputSchema["properties"]["repIndicators"]["type"] = "array";
    inputSchema["properties"]["repIndicators"]["description"] = "Repressor indicators";
    
    Json::Value required(Json::arrayValue);
    required.append("filename");
    required.append("coopMat");
    required.append("repIndicators");
    inputSchema["required"] = required;
    
    schema["inputSchema"] = inputSchema;
    return schema;
}

Json::Value ExprParLoadTool::execute(const Json::Value& parameters) {
    try {
        if (!parameters.isMember("filename") || !parameters.isMember("coopMat") ||
            !parameters.isMember("repIndicators")) {
            throw std::runtime_error("Missing required fields: filename, coopMat, repIndicators");
        }
        
        std::string filename = parameters["filename"].asString();
        IntMatrix coopMat = Serialization::jsonToIntMatrix(parameters["coopMat"]);
        
        // Convert boolean array from JSON
        std::vector<bool> repIndicators;
        for (const auto& val : parameters["repIndicators"]) {
            repIndicators.push_back(val.asBool());
        }
        
        // For demo purposes, create a mock ExprPar object since we don't have a real file
        // In a real implementation, this would use ExprPar::load()
        ExprPar par;
        
        // Mock implementation - in practice would load from file:
        // int result = par.load(filename, coopMat, repIndicators);
        // if (result != 0) throw std::runtime_error("Failed to load parameters from file");
        
        Json::Value result;
        result["success"] = true;
        result["parameters"] = Serialization::exprParToJson(par);
        result["message"] = "Parameters loaded successfully (mock implementation for file: " + filename + ")";
        
        return result;
        
    } catch (const std::exception& e) {
        Json::Value error;
        error["success"] = false;
        error["error"] = e.what();
        return error;
    }
}

ExprFuncPredictExprTool::ExprFuncPredictExprTool()
    : Tool("expr_func_predict_expr", "Predict expression using ExprFunc") {}

Json::Value ExprFuncPredictExprTool::getSchema() const {
    Json::Value schema;
    schema["name"] = getName();
    schema["description"] = getDescription();
    return schema;
}

Json::Value ExprFuncPredictExprTool::execute(const Json::Value& parameters) {
    Json::Value result;
    result["success"] = false;
    result["error"] = "Not implemented yet";
    return result;
}

// Initialize all tools
void initializeTools() {
    ToolRegistry& registry = ToolRegistry::getInstance();
    
    registry.registerTool(std::make_unique<ExprPredictorObjFuncTool>());
    registry.registerTool(std::make_unique<ExprParGetFreeParsTool>());
    registry.registerTool(std::make_unique<ExprPredictorTrainTool>());
    registry.registerTool(std::make_unique<ExprPredictorPredictTool>());
    registry.registerTool(std::make_unique<ExprParLoadTool>());
    registry.registerTool(std::make_unique<ExprFuncPredictExprTool>());
}

// MCPServer implementation
MCPServer::MCPServer() : running_(false), registry_(ToolRegistry::getInstance()) {}

Json::Value MCPServer::handleRequest(const Json::Value& request) {
    std::string method = request["method"].asString();
    
    if (method == "tools/list") {
        return handleListTools();
    } else if (method == "tools/get_schema") {
        std::string toolName = request["params"]["name"].asString();
        return handleGetToolSchema(toolName);
    } else if (method == "tools/call") {
        std::string toolName = request["params"]["name"].asString();
        Json::Value parameters = request["params"]["arguments"];
        return handleExecuteTool(toolName, parameters);
    } else {
        Json::Value error;
        error["error"] = "Unknown method: " + method;
        return error;
    }
}

Json::Value MCPServer::handleListTools() {
    Json::Value response;
    Json::Value tools(Json::arrayValue);
    
    for (const std::string& name : registry_.getToolNames()) {
        Tool* tool = registry_.getTool(name);
        Json::Value toolInfo;
        toolInfo["name"] = name;
        toolInfo["description"] = tool->getDescription();
        tools.append(toolInfo);
    }
    
    response["tools"] = tools;
    return response;
}

Json::Value MCPServer::handleGetToolSchema(const std::string& toolName) {
    Tool* tool = registry_.getTool(toolName);
    if (!tool) {
        Json::Value error;
        error["error"] = "Tool not found: " + toolName;
        return error;
    }
    
    return tool->getSchema();
}

Json::Value MCPServer::handleExecuteTool(const std::string& toolName, const Json::Value& parameters) {
    Tool* tool = registry_.getTool(toolName);
    if (!tool) {
        Json::Value error;
        error["error"] = "Tool not found: " + toolName;
        return error;
    }
    
    return tool->execute(parameters);
}

void MCPServer::start(int port) {
    running_ = true;
    std::cout << "MCP Server started on port " << port << std::endl;
    // For now, this is just a placeholder - a real implementation would start HTTP/WebSocket server
}

void MCPServer::stop() {
    running_ = false;
    std::cout << "MCP Server stopped" << std::endl;
}

} // namespace MCP