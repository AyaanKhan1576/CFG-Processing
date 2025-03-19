#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <iomanip>  // For setw

using namespace std;

//------------------------------------------------------------------------------
// Grammar structure: maps non-terminals to a list of productions,
// where each production is represented as a vector of symbols.
struct Grammar {
    map<string, vector<vector<string> > > productions;
    set<string> nonTerminals;
    set<string> terminals;
};

// Utility function: trim leading and trailing whitespace.
string trim(const string &str) {
    int start = str.find_first_not_of(" \t");
    if (start == -1) return "";
    int end = str.find_last_not_of(" \t");
    return str.substr(start, end - start + 1);
}

// Utility function: split a string by a given delimiter character.
vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        token = trim(token);
        if (token != "")
            tokens.push_back(token);
    }
    return tokens;
}

// Function to read the grammar from a file.
// Expected format per line: NonTerminal -> production1 | production2 | ...
Grammar readGrammar(const string &filename) {
    Grammar grammar;
    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }
    string line;
    while (getline(infile, line)) {
        line = trim(line);
        if (line == "")
            continue;
        int pos = line.find("->");
        if (pos == -1)
            continue;
        string lhs = trim(line.substr(0, pos));
        string rhs = trim(line.substr(pos + 2));
        grammar.nonTerminals.insert(lhs);
        vector<string> prodParts = split(rhs, '|');
        for (int i = 0; i < (int)prodParts.size(); i++) {
            vector<string> symbols = split(prodParts[i], ' ');
            grammar.productions[lhs].push_back(symbols);
        }
    }
    // Determine terminals: any symbol not in the non-terminals (and not epsilon).
    for (map<string, vector<vector<string> > >::iterator it = grammar.productions.begin(); 
         it != grammar.productions.end(); ++it) {
        for (int i = 0; i < (int)it->second.size(); i++) {
            for (int j = 0; j < (int)it->second[i].size(); j++) {
                string sym = it->second[i][j];
                if (grammar.nonTerminals.find(sym) == grammar.nonTerminals.end() && sym != "ε")
                    grammar.terminals.insert(sym);
            }
        }
    }
    return grammar;
}

// Print the grammar productions to the provided output stream.
void printGrammar(const Grammar &grammar, ostream &out) {
    out << "Grammar productions:\n";
    for (map<string, vector<vector<string> > >::const_iterator it = grammar.productions.begin(); 
         it != grammar.productions.end(); ++it) {
        out << it->first << " -> ";
        bool firstProd = true;
        vector<vector<string> > prods = it->second;
        for (int i = 0; i < (int)prods.size(); i++) {
            if (!firstProd)
                out << " | ";
            vector<string> production = prods[i];
            for (int j = 0; j < (int)production.size(); j++) {
                out << production[j] << " ";
            }
            firstProd = false;
        }
        out << "\n";
    }
    out << "\n";
}

// Utility: Compute the longest common prefix of a group of productions.
vector<string> longestCommonPrefix(const vector<vector<string> > &prods) {
    vector<string> prefix;
    if (prods.size() == 0)
        return prefix;
    prefix = prods[0];
    for (int i = 1; i < (int)prods.size(); i++) {
        int j = 0;
        int limit = (int)prefix.size();
        if (prods[i].size() < (unsigned)limit) {
            limit = (int)prods[i].size();
        }
        while (j < limit && prefix[j] == prods[i][j])
            j++;
        prefix.resize(j);
    }
    return prefix;
}

//------------------------------------------------------------------------------
// Step 1: Left Factoring
void leftFactoring(Grammar &grammar) {
    bool changes = true;
    int newNTCount = 1;
    while (changes) {
        changes = false;
        map<string, vector<vector<string> > > newProductions;
        vector< pair<string, vector<vector<string> > > > additions;
        
        // Iterate over each non-terminal production.
        for (map<string, vector<vector<string> > >::iterator it = grammar.productions.begin();
             it != grammar.productions.end(); ++it) {
            string nt = it->first;
            vector<vector<string> > prods = it->second;
            // Group productions by their first symbol.
            map<string, vector<vector<string> > > groups;
            for (int i = 0; i < (int)prods.size(); i++) {
                if (prods[i].size() > 0)
                    groups[prods[i][0]].push_back(prods[i]);
            }
            vector<vector<string> > newProdList;
            for (map<string, vector<vector<string> > >::iterator gt = groups.begin(); 
                 gt != groups.end(); ++gt) {
                if (gt->second.size() > 1) {
                    vector<string> lcp = longestCommonPrefix(gt->second);
                    if (lcp.size() == 0) {
                        for (int i = 0; i < (int)gt->second.size(); i++) {
                            newProdList.push_back(gt->second[i]);
                        }
                    } else {
                        changes = true;
                        // Create a new non-terminal name (e.g., A').
                        string newNT = nt + "'";
                        while (grammar.nonTerminals.find(newNT) != grammar.nonTerminals.end())
                        {
                            newNT = nt + "'" + to_string(newNTCount);
                            newNTCount++;
                        }
                        grammar.nonTerminals.insert(newNT);
                        vector<string> newProd = lcp;
                        newProd.push_back(newNT);
                        newProdList.push_back(newProd);
                        vector<vector<string> > newNTProds;
                        for (int i = 0; i < (int)gt->second.size(); i++) {
                            vector<string> remainder;
                            int lcpSize = (int)lcp.size();
                            for (int j = lcpSize; j < (int)gt->second[i].size(); j++) {
                                remainder.push_back(gt->second[i][j]);
                            }
                            if (remainder.size() == 0)
                                remainder.push_back("ε");
                            newNTProds.push_back(remainder);
                        }
                        additions.push_back(pair<string, vector<vector<string> > >(newNT, newNTProds));
                    }
                } else {
                    newProdList.push_back(gt->second[0]);
                }
            }
            newProductions[nt] = newProdList;
        }
        // Add new non-terminals produced.
        for (int i = 0; i < (int)additions.size(); i++) {
            newProductions[additions[i].first] = additions[i].second;
        }
        grammar.productions = newProductions;
    }
}

//------------------------------------------------------------------------------
// Step 2: Left Recursion Removal
void removeLeftRecursion(Grammar &grammar) {
    map<string, vector<vector<string> > > newProductions;
    set<string> newNonTerminals = grammar.nonTerminals;
    for (map<string, vector<vector<string> > >::iterator it = grammar.productions.begin();
         it != grammar.productions.end(); ++it) {
        string A = it->first;
        vector<vector<string> > prods = it->second;
        vector<vector<string> > alpha; // Productions of the form A -> Aα.
        vector<vector<string> > beta;  // Productions of the form A -> β.
        for (int i = 0; i < (int)prods.size(); i++) {
            if (prods[i].size() > 0 && prods[i][0] == A) {
                vector<string> rem;
                for (int j = 1; j < (int)prods[i].size(); j++) {
                    rem.push_back(prods[i][j]);
                }
                alpha.push_back(rem);
            } else {
                beta.push_back(prods[i]);
            }
        }
        if (alpha.size() > 0) {
            string Aprime = A + "'";
            while (newNonTerminals.find(Aprime) != newNonTerminals.end())
                Aprime = Aprime + "'";
            newNonTerminals.insert(Aprime);
            vector<vector<string> > newAProds;
            for (int i = 0; i < (int)beta.size(); i++) {
                vector<string> prod = beta[i];
                prod.push_back(Aprime);
                newAProds.push_back(prod);
            }
            newProductions[A] = newAProds;
            vector<vector<string> > newAprimeProds;
            for (int i = 0; i < (int)alpha.size(); i++) {
                vector<string> prod = alpha[i];
                prod.push_back(Aprime);
                newAprimeProds.push_back(prod);
            }
            vector<string> epsilonProd;
            epsilonProd.push_back("ε");
            newAprimeProds.push_back(epsilonProd);
            newProductions[Aprime] = newAprimeProds;
        } else {
            newProductions[A] = prods;
        }
    }
    grammar.productions = newProductions;
    grammar.nonTerminals = newNonTerminals;
}

//------------------------------------------------------------------------------
// Step 3: Compute FIRST Sets
map<string, set<string> > computeFirstSets(const Grammar &grammar) {
    map<string, set<string> > first;
    // Initialize FIRST sets for terminals.
    for (set<string>::iterator it = grammar.terminals.begin(); it != grammar.terminals.end(); ++it)
        first[*it].insert(*it);
    // Initialize FIRST sets for non-terminals to empty.
    for (set<string>::iterator it = grammar.nonTerminals.begin(); it != grammar.nonTerminals.end(); ++it)
        first[*it] = set<string>();
    bool changed = true;
    while (changed) {
        changed = false;
        // For each production A -> α.
        for (map<string, vector<vector<string> > >::const_iterator it = grammar.productions.begin(); 
             it != grammar.productions.end(); ++it) {
            string A = it->first;
            vector<vector<string> > prods = it->second;
            for (int i = 0; i < (int)prods.size(); i++) {
                bool epsilonInAll = true;
                for (int j = 0; j < (int)prods[i].size(); j++) {
                    string sym = prods[i][j];
                    if (sym == "ε") {
                        if (first[A].find("ε") == first[A].end()) {
                            first[A].insert("ε");
                            changed = true;
                        }
                        epsilonInAll = false;
                        break;
                    }
                    int before = first[A].size();
                    set<string> firstSym = first[sym];
                    // Insert all symbols except ε.
                    for (set<string>::iterator sit = firstSym.begin(); sit != firstSym.end(); ++sit) {
                        if (*sit != "ε")
                            first[A].insert(*sit);
                    }
                    if (firstSym.find("ε") == firstSym.end()) {
                        epsilonInAll = false;
                        break;
                    }
                    int after = first[A].size();
                    if (after > before)
                        changed = true;
                }
                if (epsilonInAll) {
                    if (first[A].find("ε") == first[A].end()) {
                        first[A].insert("ε");
                        changed = true;
                    }
                }
            }
        }
    }
    return first;
}

//------------------------------------------------------------------------------
// Step 4: Compute FOLLOW Sets
map<string, set<string> > computeFollowSets(const Grammar &grammar, 
    const map<string, set<string> > &first, const string &startSymbol) {
    map<string, set<string> > follow;
    // Initialize FOLLOW sets for non-terminals to empty.
    for (set<string>::iterator it = grammar.nonTerminals.begin(); it != grammar.nonTerminals.end(); ++it)
        follow[*it] = set<string>();
    // Add end marker to the start symbol.
    follow[startSymbol].insert("$");
    bool changed = true;
    while (changed) {
        changed = false;
        // For each production A -> α.
        for (map<string, vector<vector<string> > >::const_iterator it = grammar.productions.begin();
             it != grammar.productions.end(); ++it) {
            string A = it->first;
            vector<vector<string> > prods = it->second;
            for (int i = 0; i < (int)prods.size(); i++) {
                vector<string> prod = prods[i];
                for (int j = 0; j < (int)prod.size(); j++) {
                    string B = prod[j];
                    if (grammar.nonTerminals.find(B) != grammar.nonTerminals.end()) {
                        int before = follow[B].size();
                        bool epsilonAll = true;
                        for (int k = j + 1; k < (int)prod.size(); k++) {
                            string symbol = prod[k];
                            set<string> firstSymbol = first.at(symbol);
                            for (set<string>::iterator sit = firstSymbol.begin(); sit != firstSymbol.end(); ++sit) {
                                if (*sit != "ε")
                                    follow[B].insert(*sit);
                            }
                            if (firstSymbol.find("ε") == firstSymbol.end()) {
                                epsilonAll = false;
                                break;
                            }
                        }
                        if (epsilonAll) {
                            set<string> followA = follow[A];
                            for (set<string>::iterator sit = followA.begin(); sit != followA.end(); ++sit) {
                                follow[B].insert(*sit);
                            }
                        }
                        if ((int)follow[B].size() > before)
                            changed = true;
                    }
                }
            }
        }
    }
    return follow;
}

//------------------------------------------------------------------------------
// Step 5: Construct the LL(1) Parsing Table.
// The table is represented as a mapping: non-terminal -> (terminal -> production).
map<string, map<string, vector<string> > > constructParsingTable(
    const Grammar &grammar,
    const map<string, set<string> > &first,
    const map<string, set<string> > &follow) {
    
    map<string, map<string, vector<string> > > table;
    for (map<string, vector<vector<string> > >::const_iterator it = grammar.productions.begin();
         it != grammar.productions.end(); ++it) {
        string A = it->first;
        vector<vector<string> > prods = it->second;
        for (int i = 0; i < (int)prods.size(); i++) {
            set<string> firstAlpha;
            bool epsilonInAll = true;
            for (int j = 0; j < (int)prods[i].size(); j++) {
                string sym = prods[i][j];
                if (sym == "ε") {
                    firstAlpha.insert("ε");
                    epsilonInAll = false;
                    break;
                }
                set<string> firstSym = first.at(sym);
                for (set<string>::iterator sit = firstSym.begin(); sit != firstSym.end(); ++sit) {
                    if (*sit != "ε")
                        firstAlpha.insert(*sit);
                }
                if (firstSym.find("ε") == firstSym.end()) {
                    epsilonInAll = false;
                    break;
                }
            }
            if (epsilonInAll)
                firstAlpha.insert("ε");
            // For each terminal in FIRST(α) (except ε), add production to table.
            for (set<string>::iterator fit = firstAlpha.begin(); fit != firstAlpha.end(); ++fit) {
                if (*fit != "ε")
                    table[A][*fit] = prods[i];
            }
            // If ε is in FIRST(α), add production for each terminal in FOLLOW(A).
            if (firstAlpha.find("ε") != firstAlpha.end()) {
                set<string> followA = follow.at(A);
                for (set<string>::iterator fit = followA.begin(); fit != followA.end(); ++fit) {
                    table[A][*fit] = prods[i];
                }
            }
        }
    }
    return table;
}

//------------------------------------------------------------------------------
// Helper function to print a set mapping (e.g., FIRST or FOLLOW sets)
void printSetMap(const map<string, set<string> > &sets, const string &setName, ostream &out) {
    out << setName << " Sets:\n";
    for (map<string, set<string> >::const_iterator it = sets.begin(); it != sets.end(); ++it) {
        out << it->first << " : { ";
        set<string> innerSet = it->second;
        for (set<string>::iterator sit = innerSet.begin(); sit != innerSet.end(); ++sit)
            out << *sit << " ";
        out << "}\n";
    }
    out << "\n";
}

//------------------------------------------------------------------------------
// Helper function to print the LL(1) parsing table in a formatted table layout.
void printParsingTableFormatted(const map<string, map<string, vector<string> > > &table, 
                                const set<string> &nonTerminals,
                                const set<string> &terminals,
                                ostream &out) {
    // Prepare column headers: terminals plus the end-marker "$".
    set<string> columns = terminals;
    columns.insert("$");
    
    // Compute the width for the first column ("Non-Terminal").
    int ntWidth = (int)string("Non-Terminal").size();
    for (set<string>::const_iterator it = nonTerminals.begin(); it != nonTerminals.end(); ++it) {
        if ((int)it->size() > ntWidth)
            ntWidth = it->size();
    }
    
    // Compute column widths for each terminal column.
    map<string, int> colWidth;
    for (set<string>::const_iterator it = columns.begin(); it != columns.end(); ++it)
        colWidth[*it] = (int)it->size();
    
    // Update widths based on cell contents.
    for (set<string>::const_iterator ntit = nonTerminals.begin(); ntit != nonTerminals.end(); ++ntit) {
        for (set<string>::const_iterator colit = columns.begin(); colit != columns.end(); ++colit) {
            string cell = "";
            if (table.find(*ntit) != table.end() && table.find(*ntit)->second.find(*colit) != table.find(*ntit)->second.end()) {
                vector<string> prod = table.find(*ntit)->second.find(*colit)->second;
                for (int i = 0; i < (int)prod.size(); i++) {
                    cell += prod[i] + " ";
                }
                if (cell != "")
                    cell.erase(cell.end()-1); // remove trailing space.
            }
            if ((int)cell.size() > colWidth[*colit])
                colWidth[*colit] = cell.size();
        }
    }
    
    // Print header row.
    out << "| " << setw(ntWidth) << left << "Non-Terminal" << " |";
    for (set<string>::const_iterator colit = columns.begin(); colit != columns.end(); ++colit)
        out << " " << setw(colWidth[*colit]) << left << *colit << " |";
    out << "\n";
    
    // Print separator.
    out << "|-" << string(ntWidth, '-') << "-|";
    for (set<string>::const_iterator colit = columns.begin(); colit != columns.end(); ++colit)
        out << "-" << string(colWidth[*colit], '-') << "-|";
    out << "\n";
    
    // Print each row.
    for (set<string>::const_iterator ntit = nonTerminals.begin(); ntit != nonTerminals.end(); ++ntit) {
        out << "| " << setw(ntWidth) << left << *ntit << " |";
        for (set<string>::const_iterator colit = columns.begin(); colit != columns.end(); ++colit) {
            string cell = "";
            if (table.find(*ntit) != table.end() && table.find(*ntit)->second.find(*colit) != table.find(*ntit)->second.end()) {
                vector<string> prod = table.find(*ntit)->second.find(*colit)->second;
                for (int i = 0; i < (int)prod.size(); i++) {
                    cell += prod[i] + " ";
                }
                if (cell != "")
                    cell.erase(cell.end()-1);
            }
            out << " " << setw(colWidth[*colit]) << left << cell << " |";
        }
        out << "\n";
    }
    out << "\n";
}

//------------------------------------------------------------------------------
// Main function.
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <grammar_file>\n";
        return 1;
    }
    
    // Use an ostringstream to accumulate all outputs.
    ostringstream output;
    
    // Step 0: Read the original grammar from file.
    Grammar grammar = readGrammar(argv[1]);
    output << "Original Grammar:\n";
    printGrammar(grammar, output);
    
    // Step 1: Left Factoring.
    leftFactoring(grammar);
    output << "Grammar after Left Factoring:\n";
    printGrammar(grammar, output);
    
    // Step 2: Left Recursion Removal.
    removeLeftRecursion(grammar);
    output << "Grammar after Left Recursion Removal:\n";
    printGrammar(grammar, output);
    
    // For FIRST/FOLLOW computations, assume the first non-terminal is the start symbol.
    string startSymbol = *(grammar.nonTerminals.begin());
    
    // Step 3: Compute FIRST sets.
    map<string, set<string> > first = computeFirstSets(grammar);
    printSetMap(first, "FIRST", output);
    
    // Step 4: Compute FOLLOW sets.
    map<string, set<string> > follow = computeFollowSets(grammar, first, startSymbol);
    printSetMap(follow, "FOLLOW", output);
    
    // Step 5: Construct the LL(1) Parsing Table.
    map<string, map<string, vector<string> > > table = constructParsingTable(grammar, first, follow);
    
    // Print the formatted parsing table.
    output << "LL(1) Parsing Table:\n";
    printParsingTableFormatted(table, grammar.nonTerminals, grammar.terminals, output);
    
    // Print the accumulated output to the terminal.
    cout << output.str();
    
    // Save the output to a text file.
    ofstream fout("output_log.txt");
    if (!fout) {
        cerr << "Error opening output_log.txt for writing.\n";
        return 1;
    }
    fout << output.str();
    fout.close();
    
    return 0;
}
