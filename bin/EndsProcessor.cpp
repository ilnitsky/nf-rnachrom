#include<iostream>
#include<map>
#include<string>
#include<vector>
#include<fstream>
#include<unistd.h>
#include<fcntl.h>
#include <cstring>
#include <stdio.h>

using namespace std;

vector<string> alphabet = {"A", "G", "C", "T", "N"};
vector<string> first_read = { "", "", "", "" };
vector<string> second_read = { "", "", "", "" };
vector<string> rna_read = { "", "", "", "" };
vector<string> dna_read = { "", "", "", "" };
std::vector<std::string> split_id = {};
std::vector<std::string> empty_vec = {};
int i,j,k;
string dnabuffer = "";
string rnabuffer = "";
int success;

void print1f_simple(int filei, std::vector<string> read, std::string &wbuffer) {
    if (!read.empty())
        wbuffer += read[0] + "\n" + read[1] + "\n" + read[2] + "\n" + read[3] + "\n" ;
    if (wbuffer.length() > 10000000 || wbuffer.length() != 0 && read.empty()) {
        /*if (read.empty()) {
                                        cout << "writing " << wbuffer.length() << " bytes..." << endl;
                                }*/
                                //coute << wbuffer.c_str() << " " << endl;
                                write(filei, wbuffer.c_str(), wbuffer.length());
        //int rc = write(filei, wbuffer.c_str(), wbuffer.length());
                                //char *a ="abc";
                                //int rc = write(5, "aaa\n", 4);

                                /*if (rc == -1)
                                        perror("");
                                cout << rc << " " << errno <<endl;*/
        wbuffer="";
    }
}

void substitude_rear(vector<string> &read, string check, string sub, int &success){
  //cout << read[1] << " " << read[3] << " " << check << " " << sub << endl;
        string modread = read[1].substr(read[1].length() - check.length(), check.length() );
        //cout << "modread: " << modread << "; check: " << check << endl;
        if (modread  == check) {
    read[1] = read[1].substr(0, read[1].length() - check.length() );
                //cout << read[1] << endl;
                if (sub.length() >= check.length())  {
            //cout << read[1] << endl;
                  //cout << read[3].length() - sub.length() + check.length() << " " << sub.length() - check.length() << endl;
      read[3] += read[3].substr(read[3].length() - sub.length() + check.length(), sub.length() - check.length());
                else
                  read[3] = read[3].substr(0, read[3].length() + sub.length() - check.length());
                read[1] += sub;
        }
        else success = 0;
}

void substitude_front(vector<string> &read, string check, string sub, int &success){
  if (read[1].substr(0, check.length()) == check) {
    read[1] = read[1].substr(check.length(), read[1].length());
                if (sub.length() >= check.length()) {
                  string aux;
                  aux = read[3].substr(0, sub.length() - check.length());
                        read[3] = aux + read[3];
                }
                else
                  read[3] = read[3].substr(check.length() - sub.length(), read[3].length());
                read[1] = sub + read[1];
        }
        else success = 0;
}

void add_rear(vector<string> &read, string addition) {
  read[1] += addition;
        read[3] += read[3].substr(read[3].length() - addition.length(), addition.length());
}

void add_front(vector<string> &read, string addition) {
  read[1] = addition + read[1];
        string aux;
        aux = read[3].substr(0, addition.length());
        read[3] = aux + read[3];
}

void remove_rear(vector<string> &read, int n, int &success) {
        if (read[1].length() >= n)  {
    read[1] = read[1].substr(read[1].length() - n, n);
          read[3] = read[3].substr(read[3].length() - n, n);
        }
        else success = 0;
}

void remove_front(vector<string> &read, int n, int &success) {
  //cout << read[1].substr(n+1, read[1].length()) << endl;
        //cout << read[3].substr(n, read[3].length()) << endl;
        if (read[1].length() >= n)  {
    read[1] = read[1].substr(n, read[1].length());
    read[3] = read[3].substr(n, read[3].length());
        }
        else success = 0;
        //cout << read[3] << endl;
}

// function for spliting a line on a character
size_t split(const std::string &txt, std::vector<std::string> &strs, char ch)  {
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while( pos != std::string::npos ) {
        strs.push_back( txt.substr( initialPos, pos - initialPos ) );
        initialPos = pos + 1;

        pos = txt.find( ch, initialPos );
    }

    // Add the last one
    strs.push_back( txt.substr( initialPos, std::min( pos, txt.size() ) - initialPos + 1 ) );

    return strs.size();
}

// parser for a description sequence, changes the value of success
void parser(string describe_seq,  bool DNA_part) {
    // defining iterator for description sequnce, logical variables for states, and auxillary stings for operators
    int i;
                int deletion_number;
    bool front = true;
    bool rear = false;
    bool addition = false;
    bool deletion = false;
    bool substitution = false;
                bool sub_seq_start = false;
                string deletion_string;
    string addition_string = "";
    string check_string = "";
    string sub_string = "";
                //cout << "parsing starts" << endl;
    //  parsing description sequence
    for (i=0; i<describe_seq.length(); i++) {
                                //cout << "i = " << i << endl;
       /* string allseq = restofline;
        string allq = restoflineq;*/
        // states
        if (describe_seq[i] == '*' || describe_seq[i] == '.') {
                                                 front = false;
                                                 rear = true;
        }
                                if (describe_seq[i] == '+') addition = true;
                                if (describe_seq[i] == '-') deletion = true;
                                if (describe_seq[i] == 's') substitution = true;
                                // operators
        if (addition && describe_seq[i] != '+' && describe_seq[i] != '[' && describe_seq[i] != ']') addition_string += describe_seq[i];
                                if (deletion && describe_seq[i] != '-' && describe_seq[i] != '[' && describe_seq[i] != ']') deletion_string += describe_seq[i];
                                if (substitution && describe_seq[i] != 's' && describe_seq[i] != '[' && describe_seq[i] != '|' && describe_seq[i] != ']' && !sub_seq_start) check_string += describe_seq[i];
                                if (substitution && describe_seq[i] == '|') sub_seq_start = true;
                                if (substitution && describe_seq[i] != 's' && describe_seq[i] != '[' && describe_seq[i] != '|' && describe_seq[i] != ']' && sub_seq_start) sub_string += describe_seq[i];
                                if (describe_seq[i] == ']') {
                                                //cout << "] met" << endl;
                                                if (addition)  {
                                                        if (front) {
                                                          if (DNA_part) add_front(first_read, addition_string);
                                                                else add_front(second_read, addition_string);
                                                        }
                                                        else {
                                                          if (DNA_part) add_rear(first_read, addition_string);
                                                                else add_rear(second_read, addition_string);
                                                        }
                                                }
                                                if (deletion)  {
                                                  deletion_number = stoi(deletion_string);
                                                        if (front) {
                                                          if (DNA_part) remove_front(first_read, deletion_number, success);
                                                                else remove_front(second_read, deletion_number, success);
                                                        }
                                                        else {
                                                          if (DNA_part) remove_rear(first_read, deletion_number, success);
                                                                else remove_rear(second_read, deletion_number, success);
                                                        }
                                                }
                                                if (substitution)  {
                                                  //cout << "subs begun; front value = " << front << " DNA_part value = " << DNA_part << endl;
                                                  if (front)  {
                                                          if (DNA_part) substitude_front(first_read, check_string, sub_string, success);
                                                                else substitude_front(first_read, check_string, sub_string, success);
                                                        }
                                                        else  {
                                                          //cout << "check string = " << check_string << " sub_string = " << sub_string << endl;
                                                                //cout << "first read seq " << first_read[1] << " first read qual " << first_read[3] << endl;
                                                          if (DNA_part) substitude_rear(first_read, check_string, sub_string, success);
                                                                else substitude_rear(second_read, check_string, sub_string, success);
                                                        }
                                                }
                                                addition = false;
                                                deletion = false;
                                                substitution = false;
                                                sub_seq_start = false;
                                        }


    }
}

int main(int argc, char *argv[]) {
    map<string, int>count_dictd2r;
                map<string, int>count_dictd2f;
                map<string, int>count_dictd3f;
                map<string, int>count_dictd3r;
                map<string, int>count_dictr2r;
    map<string, int>count_dictr2f;
    map<string, int>count_dictr3f;
    map<string, int>count_dictr3r;
                string oligos_out;
                //for(i=0; i<argc; i++) cout << argv[i] << endl;
                string dna_check_front = "GC";
                string dna_sub_front = "";
                string dna_check_rear = "AG";
                string dna_sub_rear = "AGCT";
                string rna_add_front = "AA";
                string rna_check_rear = "TT";
                string rna_sub_rear = "";
                //cout << "Generating dicts" << endl;
    //Adding the elements
                for (i=0; i< 5; i++)
                                for(j=0; j<5; j++){
                                                count_dictd2f[alphabet[i] + alphabet[j]] = 0;
                                          count_dictd2r[alphabet[i] + alphabet[j]] = 0;
                                                count_dictr2f[alphabet[i] + alphabet[j]] = 0;
            count_dictr2r[alphabet[i] + alphabet[j]] = 0;
                                }

                for (i=0; i< 5; i++)
        for(j=0; j<5; j++)
                                  for(k=0; k<5; k++){
                                          count_dictd3f[alphabet[i] + alphabet[j] + alphabet[k]] = 0;
                                                count_dictd3r[alphabet[i] + alphabet[j] + alphabet[k]] = 0;
                                                count_dictr3f[alphabet[i] + alphabet[j] + alphabet[k]] = 0;
            count_dictr3r[alphabet[i] + alphabet[j] + alphabet[k]] = 0;
                                  }
                //cout << "Dicts are initiated" << endl;
                //cout << alphabet[i] + alphabet[j] << endl;
                //cout << count_dict["AA"] << endl;
                //Traversing through the map elements
    /*for (auto element :count_dict){

        //element.first represents key
        cout<<element.first<<" is the capital of ";

        //element.second represents value
        cout<<element.second<<endl;
    }
    */
    // Removing an element

    // Size of Map
    //cout<<"The size of Map is :"<< count_dict.size();
                string inputfilepath1 = argv[1];
                string inputfilepath2 = argv[2];
                string desc_seq_dna = argv[3];
                string desc_seq_rna = argv[4];
                cout << "End sequence processor has started" << endl;
                //cout << "desc_seq DNA = " << desc_seq_dna <<  "desc_seq RNA = " << desc_seq_rna << endl;
                string s1 = inputfilepath1;
                s1.append("_RS");
                string s2 = inputfilepath2;
                s2.append("_RS");
                string s3 = inputfilepath1;
                s3.append("_last_oligos.tsv");
                fstream infile1;
                fstream infile2;
                infile1.open(inputfilepath1, std::ios::in);
                infile2.open(inputfilepath2, std::ios::in);
//              const int dnafile = open(s1.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);
//              const int rnafile = open(s2.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);
                const int dnafile = open(s1.c_str(), O_CREAT | O_WRONLY, 0644);
                const int rnafile = open(s2.c_str(), O_CREAT | O_WRONLY, 0644);
                const int oligos = open(s3.c_str(), O_CREAT | O_WRONLY, 0644);
                if (access (s1.c_str(), W_OK) ||  access (s2.c_str(), W_OK) ||  access (s2.c_str(), W_OK)) {
        cerr << "Cannot write to output file" << endl;
        exit(1);
    }
                string line1, line2;
                int line_cnt_mod = 0;
                while(getline(infile1, line1) && getline(infile2, line2)) {
                                //cout << line_cnt_mod << endl;
                                first_read[line_cnt_mod] = line1;
        second_read[line_cnt_mod] = line2;
                                if (line_cnt_mod == 3) {
                                                split( first_read[0], split_id, ' ' );
                                                first_read[0] = split_id[0];
                                                split_id.clear();
                                                split( second_read[0], split_id, ' ' );
            second_read[0] = split_id[0];
            split_id.clear();
                                                //icout << first_read[0] << " " << second_read[0] << endl;
                                                /*
                                                if (first_read[0] != second_read[0]) {
                    cerr << "input files are not sorted, different read IDs in the lines" << endl;
                    exit(1);
            }*/
                                                //cout << "before counting iteration" << endl;
                                                count_dictd2f[first_read[1].substr(0, 2)] += 1;
                                                count_dictd2r[first_read[1].substr(first_read[1].length() - 2, 2)] += 1;
                                                count_dictr2f[second_read[1].substr(0, 2)] += 1;
                                                count_dictr2r[second_read[1].substr(second_read[1].length() - 2, 2)] += 1;
                                                count_dictd3f[first_read[1].substr(0, 3)] += 1;
            count_dictd3r[first_read[1].substr(first_read[1].length() - 3, 3)] += 1;
            count_dictr3f[second_read[1].substr(0, 3)] += 1;
            count_dictr3r[second_read[1].substr(second_read[1].length() - 3, 3)] += 1;
                                                //cout << "after counting iteration/before ops" << endl;
                                                //substitude_rear(first_read, dna_check_rear, dna_sub_rear);
                                                //cout << "after sub rear" << endl;
                                                //add_front(second_read, rna_add_front);
                                                //remove_front(second_read, 1);
                                                //cout << second_read[3] << endl;
                                                //cout << "after add front" << endl;
                                                //substitude_rear(second_read, rna_check_rear, rna_sub_rear);
                                                //substitude_front(first_read, dna_check_front, dna_sub_front);
                                                //cout << "ops finished" << endl;

                                                /*if (first_read[1].substr(first_read[1].length() - 2, 2) == "AG") {
                                                                first_read[1].append("CT");
                                                                first_read[3].append(first_read[3].substr(first_read[3].length() - 2, 2));
                                                }*/

                                                /*for (i=0; i<first_read.size(); i++)
                                                                cout << first_read[i] << endl;
                                                for (i=0; i<second_read.size(); i++)
                cout << second_read[i] << endl;*/
                                                success = 1;
                                                //cout << "before dna parser" << endl;
                                                parser(desc_seq_dna, true);
                                                parser(desc_seq_rna, false);
                                                if (success)  {
                                                  print1f_simple(dnafile, first_read, dnabuffer);
                                                  print1f_simple(rnafile, second_read, rnabuffer);
                                                }
                                                first_read = { "", "", "", "" };
                                                second_read = { "", "", "", "" };
                                                line_cnt_mod = 0;
                                }
                                line_cnt_mod += 1;
                }
                //cout << empty_vec.empty() << " " << dnabuffer.length() << endl;
                line_cnt_mod = 0;
                print1f_simple(dnafile, empty_vec, dnabuffer);
    print1f_simple(rnafile, empty_vec, rnabuffer);
                close(dnafile);
                close(rnafile);
                for (auto element :count_dictd2r)  {
                                oligos_out = element.first + "\tDNA_3'end_last_2:\t" + to_string(element.second) + "\n";
                                write(oligos, oligos_out.c_str(), oligos_out.length());
                }
                for (auto element :count_dictd2f)  {
                                oligos_out = element.first + "\tDNA_5'end_last_2:\t" + to_string(element.second) + "\n";
                                write(oligos, oligos_out.c_str(), oligos_out.length());
                }
                for (auto element :count_dictr2r)  {
                                oligos_out = element.first + "\tRNA_3'end_last_2:\t" + to_string(element.second) + "\n";
                                write(oligos, oligos_out.c_str(), oligos_out.length());
                }
                for (auto element :count_dictr2f)  {
                                oligos_out = element.first + "\tRNA_5'end_last_2:\t" + to_string(element.second) + "\n";
                                write(oligos, oligos_out.c_str(), oligos_out.length());
                }
                for (auto element :count_dictd3r)  {
                                oligos_out = element.first + "\tDNA_3'end_last_3:\t" + to_string(element.second) + "\n";
                                write(oligos, oligos_out.c_str(), oligos_out.length());
                }
    for (auto element :count_dictd3f)  {
                                oligos_out = element.first + "\tDNA_5'end_last_3:\t" + to_string(element.second) + "\n";
        write(oligos, oligos_out.c_str(), oligos_out.length());
                }
    for (auto element :count_dictr3r)  {
                                oligos_out = element.first + "\tRNA_3'end_last_3:\t" + to_string(element.second) + "\n";
        write(oligos, oligos_out.c_str(), oligos_out.length());
                }
    for (auto element :count_dictr3f)  {
                                oligos_out = element.first + "\tRNA_5'end_last_3:\t" + to_string(element.second) + "\n";
        write(oligos, oligos_out.c_str(), oligos_out.length());
                }
                close(oligos);
    // Clear the Map
    count_dictd2r.clear();
                count_dictd2f.clear();
                count_dictr2r.clear();
    count_dictr2f.clear();
                count_dictd3r.clear();
    count_dictd3f.clear();
    count_dictr3r.clear();
    count_dictr3f.clear();
                return 0;
}
