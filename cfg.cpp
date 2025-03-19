#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

using namespace std;

// Structure to represent a production rule
struct Production {
    string lhs;                // Left-hand side (non-terminal)
    vector<string> rhs;        // Right-hand side (sequence of terminals and non-terminals)
};

// Structure to represent a CFG
struct Grammar {
    set<string> terminals;
    set<string> nonTerminals;
    map<string, vector<vector<string>>> productions; // Non-terminal -> list of alternative productions
    string startSymbol;
};

// Class for processing CFG
class CFGProcessor {
private:
    Grammar grammar;
    map<string, set<string>> firstSets;
    map<string, set<string>> followSets;
    map<pair<string, string>, vector<string>> parseTable;
    bool isTerminal(const string& symbol);
    bool isNonTerminal(const string& symbol);
    set<string> computeFirstOfString(const vector<string>& symbols);
    ofstream outputFile;

public:
    CFGProcessor(const string& filename, const string& outputFilename);
    ~CFGProcessor();
    void displayGrammar(const Grammar& g);
    void performLeftFactoring();
    void eliminateLeftRecursion();
    void computeFirstSets();
    void computeFollowSets();
    void constructParseTable();
    void displayResults();
};

// Constructor to read grammar from file and open output file
CFGProcessor::CFGProcessor(const string& filename, const string& outputFilename) {
    // Open the output file
    outputFile.open(outputFilename);
    if (!outputFile.is_open()) {
        cerr << "Error opening output file: " << outputFilename << endl;
        exit(1);
    }

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        outputFile.close();
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Parse the production rule
        size_t arrowPos = line.find("->");
        if (arrowPos == string::npos) continue;

        string lhs = line.substr(0, arrowPos);
        string rhs = line.substr(arrowPos + 2);

        // Trim whitespace
        lhs.erase(0, lhs.find_first_not_of(" \t"));
        lhs.erase(lhs.find_last_not_of(" \t") + 1);
        rhs.erase(0, rhs.find_first_not_of(" \t"));
        rhs.erase(rhs.find_last_not_of(" \t") + 1);

        // Add LHS to non-terminals
        grammar.nonTerminals.insert(lhs);

        // If this is the first production, set as start symbol
        if (grammar.startSymbol.empty()) {
            grammar.startSymbol = lhs;
        }

        // Parse the RHS alternatives (separated by '|')
        istringstream rhsStream(rhs);
        string alternative;
        while (getline(rhsStream, alternative, '|')) {
            // Trim whitespace
            alternative.erase(0, alternative.find_first_not_of(" \t"));
            alternative.erase(alternative.find_last_not_of(" \t") + 1);

            vector<string> symbols;
            istringstream symbolStream(alternative);
            string symbol;
            while (symbolStream >> symbol) {
                symbols.push_back(symbol);
                
                // We'll determine terminals after all non-terminals are known
            }
            
            // Handle the empty alternative (epsilon)
            if (symbols.empty()) {
                symbols.push_back("epsilon");
            }
            
            grammar.productions[lhs].push_back(symbols);
        }
    }
    file.close();
    
    // Now determine terminals after all non-terminals are identified
    grammar.terminals.insert("epsilon");
    
    // Process all symbols in productions to identify terminals
    for (const auto& entry : grammar.productions) {
        for (const auto& prod : entry.second) {
            for (const auto& symbol : prod) {
                if (symbol != "epsilon" && grammar.nonTerminals.find(symbol) == grammar.nonTerminals.end()) {
                    grammar.terminals.insert(symbol);
                }
            }
        }
    }
}

// Destructor to close the output file
CFGProcessor::~CFGProcessor() {
    if (outputFile.is_open()) {
        outputFile.close();
    }
}

// Check if a symbol is a terminal
bool CFGProcessor::isTerminal(const string& symbol) {
    return grammar.terminals.find(symbol) != grammar.terminals.end();
}

// Check if a symbol is a non-terminal
bool CFGProcessor::isNonTerminal(const string& symbol) {
    return grammar.nonTerminals.find(symbol) != grammar.nonTerminals.end();
}

// Display the grammar
void CFGProcessor::displayGrammar(const Grammar& g) {
    cout << "Grammar:" << endl;
    outputFile << "Grammar:" << endl;
    for (const auto& entry : g.productions) {
        cout << entry.first << " -> ";
        outputFile << entry.first << " -> ";
        for (size_t i = 0; i < entry.second.size(); ++i) {
            if (i > 0) {
                cout << " | ";
                outputFile << " | ";
            }
            for (const auto& symbol : entry.second[i]) {
                cout << symbol << " ";
                outputFile << symbol << " ";
            }
        }
        cout << endl;
        outputFile << endl;
    }
    cout << endl;
    outputFile << endl;
}

// Perform left factoring on the grammar
void CFGProcessor::performLeftFactoring() {
    Grammar newGrammar = grammar;
    bool factored = false;
    
    do {
        factored = false;
        Grammar tempGrammar = newGrammar;
        newGrammar.productions.clear();
        
        for (const auto& entry : tempGrammar.productions) {
            string nonTerminal = entry.first;
            vector<vector<string>> productions = entry.second;
            
            // Map prefix -> list of productions with that prefix
            map<string, vector<vector<string>>> prefixMap;
            for (const auto& prod : productions) {
                if (prod.empty()) continue;
                string prefix = prod[0];
                prefixMap[prefix].push_back(prod);
            }
            
            bool localFactored = false;
            for (const auto& prefixEntry : prefixMap) {
                // If a prefix appears in multiple productions, factor it out
                if (prefixEntry.second.size() > 1) {
                    localFactored = true;
                    factored = true;
                    
                    // Create a new non-terminal
                    string newNonTerminal = nonTerminal + "_" + to_string(newGrammar.nonTerminals.size());
                    newGrammar.nonTerminals.insert(newNonTerminal);
                    
                    // Add the factored production
                    vector<string> prefixProd = {prefixEntry.first, newNonTerminal};
                    newGrammar.productions[nonTerminal].push_back(prefixProd);
                    
                    // Add the new productions for the new non-terminal
                    for (const auto& prod : prefixEntry.second) {
                        vector<string> newProd;
                        // Skip the common prefix
                        for (size_t i = 1; i < prod.size(); ++i) {
                            newProd.push_back(prod[i]);
                        }
                        // If the result is empty, add epsilon
                        if (newProd.empty()) {
                            newProd.push_back("epsilon");
                        }
                        newGrammar.productions[newNonTerminal].push_back(newProd);
                    }
                } else {
                    // Keep productions without common prefixes
                    for (const auto& prod : prefixEntry.second) {
                        newGrammar.productions[nonTerminal].push_back(prod);
                    }
                }
            }
            
            // If no factoring was done for this non-terminal, keep the original productions
            if (!localFactored) {
                newGrammar.productions[nonTerminal] = productions;
            }
        }
    } while (factored);
    
    grammar = newGrammar;
    cout << "Grammar after Left Factoring:" << endl;
    outputFile << "Grammar after Left Factoring:" << endl;
    displayGrammar(grammar);
}

// Eliminate left recursion from the grammar
void CFGProcessor::eliminateLeftRecursion() {
    Grammar newGrammar = grammar;
    newGrammar.productions.clear();
    
    // Get a vector of non-terminals in order
    vector<string> nonTerminals;
    for (const auto& nt : grammar.nonTerminals) {
        nonTerminals.push_back(nt);
    }
    
    // Eliminate left recursion using the algorithm
    for (size_t i = 0; i < nonTerminals.size(); ++i) {
        string Ai = nonTerminals[i];
        
        // Substitute earlier non-terminals
        vector<vector<string>> newProductions;
        for (const auto& gamma : grammar.productions[Ai]) {
            if (!gamma.empty() && isNonTerminal(gamma[0]) && gamma[0] != Ai) {
                string Aj = gamma[0];
                if (find(nonTerminals.begin(), nonTerminals.begin() + i, Aj) != nonTerminals.begin() + i) {
                    // Substitute Aj with its productions
                    for (const auto& delta : newGrammar.productions[Aj]) {
                        vector<string> newProd = delta;
                        // Append the rest of gamma
                        for (size_t k = 1; k < gamma.size(); ++k) {
                            newProd.push_back(gamma[k]);
                        }
                        newProductions.push_back(newProd);
                    }
                } else {
                    newProductions.push_back(gamma);
                }
            } else {
                newProductions.push_back(gamma);
            }
        }
        
        // Eliminate immediate left recursion
        vector<vector<string>> alphaProductions;
        vector<vector<string>> betaProductions;
        
        for (const auto& prod : newProductions) {
            if (!prod.empty() && prod[0] == Ai) {
                // This is a left-recursive production Ai -> Ai alpha
                vector<string> alpha(prod.begin() + 1, prod.end());
                alphaProductions.push_back(alpha);
            } else {
                // This is a non-left-recursive production Ai -> beta
                betaProductions.push_back(prod);
            }
        }
        
        // If there are left-recursive productions, eliminate them
        if (!alphaProductions.empty()) {
            string newNonTerminal = Ai + "'";
            newGrammar.nonTerminals.insert(newNonTerminal);
            
            // Ai -> beta Ai'
            for (const auto& beta : betaProductions) {
                vector<string> newProd = beta;
                newProd.push_back(newNonTerminal);
                newGrammar.productions[Ai].push_back(newProd);
            }
            
            // If no beta productions, add Ai -> Ai'
            if (betaProductions.empty()) {
                newGrammar.productions[Ai].push_back({newNonTerminal});
            }
            
            // Ai' -> alpha Ai' | epsilon
            for (const auto& alpha : alphaProductions) {
                vector<string> newProd = alpha;
                newProd.push_back(newNonTerminal);
                newGrammar.productions[newNonTerminal].push_back(newProd);
            }
            
            // Add epsilon production for Ai'
            newGrammar.productions[newNonTerminal].push_back({"epsilon"});
        } else {
            // No left recursion, just copy the productions
            newGrammar.productions[Ai] = newProductions;
        }
    }
    
    grammar = newGrammar;
    cout << "Grammar after Left Recursion Elimination:" << endl;
    outputFile << "Grammar after Left Recursion Elimination:" << endl;
    displayGrammar(grammar);
}

// Compute FIRST set for a string of symbols
set<string> CFGProcessor::computeFirstOfString(const vector<string>& symbols) {
    set<string> firstSet;
    
    if (symbols.empty()) {
        firstSet.insert("epsilon");
        return firstSet;
    }
    
    string firstSymbol = symbols[0];
    
    // If the first symbol is a terminal, FIRST is just that terminal
    if (isTerminal(firstSymbol)) {
        if (firstSymbol != "epsilon") {
            firstSet.insert(firstSymbol);
        } else if (symbols.size() > 1) {
            // If first symbol is epsilon, compute FIRST of the rest
            vector<string> restSymbols(symbols.begin() + 1, symbols.end());
            set<string> restFirst = computeFirstOfString(restSymbols);
            firstSet.insert(restFirst.begin(), restFirst.end());
        } else {
            firstSet.insert("epsilon");
        }
        return firstSet;
    }
    
    // If the first symbol is a non-terminal
    // Add all elements from FIRST(firstSymbol) except epsilon
    for (const auto& term : firstSets[firstSymbol]) {
        if (term != "epsilon") {
            firstSet.insert(term);
        }
    }
    
    // Check if FIRST(firstSymbol) contains epsilon
    bool epsilonInFirst = firstSets[firstSymbol].find("epsilon") != firstSets[firstSymbol].end();
    
    // If epsilon is in FIRST(firstSymbol) and there are more symbols, compute FIRST of the rest
    if (epsilonInFirst && symbols.size() > 1) {
        vector<string> restSymbols(symbols.begin() + 1, symbols.end());
        set<string> restFirst = computeFirstOfString(restSymbols);
        firstSet.insert(restFirst.begin(), restFirst.end());
    } else if (epsilonInFirst && symbols.size() == 1) {
        // If epsilon is in FIRST(firstSymbol) and there are no more symbols
        firstSet.insert("epsilon");
    }
    
    return firstSet;
}

// Compute FIRST sets for all non-terminals
void CFGProcessor::computeFirstSets() {
    // Initialize all FIRST sets as empty
    for (const auto& nt : grammar.nonTerminals) {
        firstSets[nt] = set<string>();
    }
    
    // For terminals, FIRST is just the terminal itself
    for (const auto& t : grammar.terminals) {
        firstSets[t] = {t};
    }
    
    bool changed;
    do {
        changed = false;
        
        for (const auto& entry : grammar.productions) {
            string nonTerminal = entry.first;
            
            for (const auto& production : entry.second) {
                set<string> productionFirst;
                
                if (production.empty() || production[0] == "epsilon") {
                    // Epsilon production
                    productionFirst.insert("epsilon");
                } else {
                    // Compute FIRST of the production
                    set<string> first = computeFirstOfString(production);
                    productionFirst.insert(first.begin(), first.end());
                }
                
                // Add to FIRST(nonTerminal)
                size_t beforeSize = firstSets[nonTerminal].size();
                firstSets[nonTerminal].insert(productionFirst.begin(), productionFirst.end());
                if (firstSets[nonTerminal].size() > beforeSize) {
                    changed = true;
                }
            }
        }
    } while (changed);
    
    // Display FIRST sets
    cout << "FIRST Sets:" << endl;
    outputFile << "FIRST Sets:" << endl;
    for (const auto& entry : firstSets) {
        if (isNonTerminal(entry.first)) {
            cout << "FIRST(" << entry.first << ") = { ";
            outputFile << "FIRST(" << entry.first << ") = { ";
            bool first = true;
            for (const auto& symbol : entry.second) {
                if (!first) {
                    cout << ", ";
                    outputFile << ", ";
                }
                cout << symbol;
                outputFile << symbol;
                first = false;
            }
            cout << " }" << endl;
            outputFile << " }" << endl;
        }
    }
    cout << endl;
    outputFile << endl;
}

// Compute FOLLOW sets for all non-terminals
void CFGProcessor::computeFollowSets() {
    // Initialize all FOLLOW sets as empty
    for (const auto& nt : grammar.nonTerminals) {
        followSets[nt] = set<string>();
    }
    
    // Remove the line that adds $ to FOLLOW(start symbol)
    // followSets[grammar.startSymbol].insert("$");
    
    bool changed;
    do {
        changed = false;
        
        for (const auto& entry : grammar.productions) {
            string nonTerminal = entry.first;
            
            for (const auto& production : entry.second) {
                for (size_t i = 0; i < production.size(); ++i) {
                    // Only interested in non-terminals in the RHS
                    if (!isNonTerminal(production[i])) continue;
                    
                    string B = production[i];
                    bool isLast = (i == production.size() - 1);
                    
                    if (isLast) {
                        // If B is the last symbol, add FOLLOW(A) to FOLLOW(B)
                        size_t beforeSize = followSets[B].size();
                        followSets[B].insert(followSets[nonTerminal].begin(), followSets[nonTerminal].end());
                        if (followSets[B].size() > beforeSize) {
                            changed = true;
                        }
                    } else {
                        // Compute FIRST of the rest of the production
                        vector<string> beta(production.begin() + i + 1, production.end());
                        set<string> firstBeta = computeFirstOfString(beta);
                        
                        // Add FIRST(beta) - {epsilon} to FOLLOW(B)
                        for (const auto& symbol : firstBeta) {
                            if (symbol != "epsilon") {
                                if (followSets[B].find(symbol) == followSets[B].end()) {
                                    followSets[B].insert(symbol);
                                    changed = true;
                                }
                            }
                        }
                        
                        // If epsilon is in FIRST(beta), add FOLLOW(A) to FOLLOW(B)
                        if (firstBeta.find("epsilon") != firstBeta.end()) {
                            size_t beforeSize = followSets[B].size();
                            followSets[B].insert(followSets[nonTerminal].begin(), followSets[nonTerminal].end());
                            if (followSets[B].size() > beforeSize) {
                                changed = true;
                            }
                        }
                    }
                }
            }
        }
    } while (changed);
    
    // Display FOLLOW sets
    cout << "FOLLOW Sets:" << endl;
    outputFile << "FOLLOW Sets:" << endl;
    for (const auto& entry : followSets) {
        cout << "FOLLOW(" << entry.first << ") = { ";
        outputFile << "FOLLOW(" << entry.first << ") = { ";
        bool first = true;
        for (const auto& symbol : entry.second) {
            if (!first) {
                cout << ", ";
                outputFile << ", ";
            }
            cout << symbol;
            outputFile << symbol;
            first = false;
        }
        cout << " }" << endl;
        outputFile << " }" << endl;
    }
    cout << endl;
    outputFile << endl;
}

// Construct the LL(1) parsing table and display in the requested format
void CFGProcessor::constructParseTable() {
    parseTable.clear();
    
    // For each production A -> α
    for (const auto& entry : grammar.productions) {
        string nonTerminal = entry.first;
        
        for (const auto& production : entry.second) {
            // Compute FIRST(α)
            set<string> firstAlpha = computeFirstOfString(production);
            
            // For each terminal 'a' in FIRST(α), add A -> α to M[A, a]
            for (const auto& terminal : firstAlpha) {
                if (terminal != "epsilon") {
                    parseTable[{nonTerminal, terminal}] = production;
                }
            }
            
            // If epsilon is in FIRST(α), for each terminal 'b' in FOLLOW(A), add A -> α to M[A, b]
            if (firstAlpha.find("epsilon") != firstAlpha.end()) {
                for (const auto& terminal : followSets[nonTerminal]) {
                    parseTable[{nonTerminal, terminal}] = production;
                }
            }
        }
    }
    
    // Format and display the parsing table as in the image
    cout << "LL(1) Parsing Table:" << endl;
    outputFile << "LL(1) Parsing Table:" << endl;
    
    // Get all terminals (excluding epsilon) for columns
    set<string> tableTerminals;
    for (const auto& term : grammar.terminals) {
        if (term != "epsilon") tableTerminals.insert(term);
    }
    // Remove the line that adds $ to tableTerminals
    // tableTerminals.insert("$"); // Add end marker
    
    // Print the table header with terminals as columns
    const int colWidth = 15;
    
    // Print top border
    cout << "+";
    outputFile << "+";
    cout << string(colWidth, '-') << "+";
    outputFile << string(colWidth, '-') << "+";
    for (const auto& term : tableTerminals) {
        cout << string(colWidth, '-') << "+";
        outputFile << string(colWidth, '-') << "+";
    }
    cout << endl;
    outputFile << endl;
    
    // Print terminal headers
    cout << "|" << setw(colWidth) << "  " << "|";
    outputFile << "|" << setw(colWidth) << "  " << "|";
    for (const auto& term : tableTerminals) {
        cout << setw(colWidth) << term << "|";
        outputFile << setw(colWidth) << term << "|";
    }
    cout << endl;
    outputFile << endl;
    
    // Print border
    cout << "+";
    outputFile << "+";
    cout << string(colWidth, '-') << "+";
    outputFile << string(colWidth, '-') << "+";
    for (const auto& term : tableTerminals) {
        cout << string(colWidth, '-') << "+";
        outputFile << string(colWidth, '-') << "+";
    }
    cout << endl;
    outputFile << endl;
    
    // Print rows for each non-terminal
    for (const auto& nt : grammar.nonTerminals) {
        cout << "|" << setw(colWidth) << nt << "|";
        outputFile << "|" << setw(colWidth) << nt << "|";
        
        for (const auto& term : tableTerminals) {
            string cellContent = "";
            if (parseTable.find({nt, term}) != parseTable.end()) {
                cellContent = nt + " -> ";
                for (const auto& symbol : parseTable[{nt, term}]) {
                    cellContent += symbol + " ";
                }
            }
            cout << setw(colWidth) << cellContent << "|";
            outputFile << setw(colWidth) << cellContent << "|";
        }
        cout << endl;
        outputFile << endl;
        
        // Print border after each row
        cout << "+";
        outputFile << "+";
        cout << string(colWidth, '-') << "+";
        outputFile << string(colWidth, '-') << "+";
        for (const auto& term : tableTerminals) {
            cout << string(colWidth, '-') << "+";
            outputFile << string(colWidth, '-') << "+";
        }
        cout << endl;
        outputFile << endl;
    }
}
// Display all results of the CFG processing
void CFGProcessor::displayResults() {
    cout << "Original Grammar:" << endl;
    outputFile << "Original Grammar:" << endl;
    displayGrammar(grammar);
    
    performLeftFactoring();
    eliminateLeftRecursion();
    computeFirstSets();
    computeFollowSets();
    constructParseTable();
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " grammar.txt output.txt" << endl;
        return 1;
    }
    
    CFGProcessor processor(argv[1], argv[2]);
    processor.displayResults();
    
    cout << "Results have been saved to " << argv[2] << endl;
    
    return 0;
}