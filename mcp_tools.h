#ifndef MCP_TOOLS_H
#define MCP_TOOLS_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <jsoncpp/json/json.h>
#include "ExprPredictor.h"

/**
 * MCP (Model Context Protocol) Tools for ExprPredictor
 * 
 * This header defines the MCP tool interface for exposing ExprPredictor
 * functions as remotely callable tools through the MCP protocol.
 */

namespace MCP {

// Forward declarations
class ToolRegistry;
class Tool;

/**
 * Base class for MCP tool implementations
 */
class Tool {
public:
    Tool(const std::string& name, const std::string& description)
        : name_(name), description_(description) {}
    
    virtual ~Tool() = default;
    
    virtual Json::Value getSchema() const = 0;
    virtual Json::Value execute(const Json::Value& parameters) = 0;
    
    const std::string& getName() const { return name_; }
    const std::string& getDescription() const { return description_; }

protected:
    std::string name_;
    std::string description_;
};

/**
 * Tool registry for managing MCP tools
 */
class ToolRegistry {
public:
    static ToolRegistry& getInstance();
    
    void registerTool(std::unique_ptr<Tool> tool);
    Tool* getTool(const std::string& name);
    std::vector<std::string> getToolNames() const;
    Json::Value getToolsSchema() const;

private:
    std::map<std::string, std::unique_ptr<Tool>> tools_;
};

/**
 * Utility functions for JSON serialization
 */
namespace Serialization {
    Json::Value vectorToJson(const std::vector<double>& vec);
    Json::Value matrixToJson(const Matrix& matrix);
    Json::Value intMatrixToJson(const IntMatrix& matrix);
    
    std::vector<double> jsonToVector(const Json::Value& json);
    Matrix jsonToMatrix(const Json::Value& json);
    IntMatrix jsonToIntMatrix(const Json::Value& json);
    
    Json::Value exprParToJson(const ExprPar& par);
    ExprPar jsonToExprPar(const Json::Value& json);
}

/**
 * MCP Tools for ExprPredictor functions
 */

class ExprPredictorTrainTool : public Tool {
public:
    ExprPredictorTrainTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

class ExprPredictorPredictTool : public Tool {
public:
    ExprPredictorPredictTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

class ExprPredictorObjFuncTool : public Tool {
public:
    ExprPredictorObjFuncTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

class ExprParGetFreeParsTool : public Tool {
public:
    ExprParGetFreeParsTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

class ExprParLoadTool : public Tool {
public:
    ExprParLoadTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

class ExprFuncPredictExprTool : public Tool {
public:
    ExprFuncPredictExprTool();
    Json::Value getSchema() const override;
    Json::Value execute(const Json::Value& parameters) override;
};

/**
 * Initialize all MCP tools
 */
void initializeTools();

/**
 * MCP Server interface
 */
class MCPServer {
public:
    MCPServer();
    
    void start(int port = 8080);
    void stop();
    
    Json::Value handleRequest(const Json::Value& request);

private:
    bool running_;
    ToolRegistry& registry_;
    
    Json::Value handleListTools();
    Json::Value handleGetToolSchema(const std::string& toolName);
    Json::Value handleExecuteTool(const std::string& toolName, const Json::Value& parameters);
};

} // namespace MCP

#endif // MCP_TOOLS_H