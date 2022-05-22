/**
 * @file finger_grammar.h
 * @author Jie Xu 
 */ 
#pragma once

#include "Common/Common.h"
#include "Graph/tree.h"
#include "Graph/node.h"

/**
 * Implement a rule in the context-free grammar.
 * 
 * Each rule contains a LHS symbol and list of RHS symbols.
 * context-free rule requires the LHS is a single nonterminal symbol, and here 
 * further require the RHS is a sequential symbol list (no branches)
 */
class FingerGrammarRule {
public:
    FingerGrammarRule() {}
    FingerGrammarRule(std::string LHS, std::vector<std::string> RHS);

    std::string _LHS;

    // Template node for each symbol in RHS, 
    // the transformation stored in node._T_p
    std::vector<std::string> _RHS; 
    std::vector<Node> _RHS_nodes;
};

/**
 * Implement a context-free finger grammar.
 * 
 * The grammar contains a list of production rules.
 */
class FingerGrammar {
public:
    FingerGrammar();
    FingerGrammar(std::string grammar_file);

    void add_rule(FingerGrammarRule rule);

    Tree get_init_graph();

    void reset();
    void undo();

    std::vector<int> get_available_rules();

    bool apply_rule(int rule_id);

    void finalize_rule_nodes();

    void print_grammar();

    std::string _name;
    
    std::vector<FingerGrammarRule> _rules;

    std::string _starting_symbol;

    std::map<std::string, Node> _template_nodes;

    Tree _graph;
    std::vector<Tree> _graph_his;

private:
    Eigen::VectorXd str_to_eigen(std::string str);
};