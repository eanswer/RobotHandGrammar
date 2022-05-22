/**
 * @file tree.h
 * @author Jie Xu 
 */ 
#pragma once

#include "Common/Common.h"
#include "Graph/node.h"

class LinearBlendSkinningObject;
class FingerGrammarRule;
class PalmGrammar;

/**
 * Implement a directed design graph (a tree). 
 * 
 * Each tree contains a list of nodes and a root node.
 */
class Tree {
public:
    Tree();
    Tree(Node* root);

    void copy_from(const Tree& tree);
    Node* copy(const Node* root);
    
    // find the a node by a rule, return the first on in dfs order
    Node* find_node(const FingerGrammarRule rule) const;
    Node* find_node(Node* node, const FingerGrammarRule rule) const;

    // replace a node by a chain, if chain has zero node, then it's equivalent
    // to deleting the node
    void replace(Node* node, std::vector<Node*> chain); 

    /**
     * get objects and handle constraints from the design graph.
     * 
     * Output:
     * objects: LinearBlendSkinningObject
     * types: node type (TERMINAL or NONTERMINAL)
     * constraints: the handle correspondences for constraints.
     */
    void get_objects_and_constraints(std::vector<LinearBlendSkinningObject*>& objects,
                    std::vector<Node::NodeType>& types,
                    std::map<std::pair<int, int>, std::pair<int, int> > &constraints);
    void get_objects_and_constraints(Node* root,
                    int parent_id,
                    Eigen::Matrix4d T, 
                    std::vector<LinearBlendSkinningObject*>& objects, 
                    std::vector<Node::NodeType>& types,
                    std::map<std::pair<int, int>, std::pair<int, int> > &constraints);

    void get_deformation_objects(std::vector<LinearBlendSkinningObject*>& objects,
                    std::vector<Node::NodeType>& types,
                    PalmGrammar* palm_grammar);
    void get_deformation_objects(Node* root,
                    int parent_id,
                    Eigen::Matrix4d T, 
                    std::vector<LinearBlendSkinningObject*>& objects, 
                    std::vector<Node::NodeType>& types,
                    PalmGrammar* palm_grammar);

    void print_graph();
    void print_graph(Node* node);
    
    Node* _root;
};