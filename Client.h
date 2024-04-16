
/*

client by ning.
This implementation is based on Feather.

*/


//*********************************************************************
#include "Server.h"

//**********************************************************************

class Client{

public:
	Client();
	Client (Server* serv, bigint *, int elem_size);
	void free_client();
	GrantComp_Info * grant_comp(CompPerm_Request* , bigint **&qq, bool);
	vector <string> find_intersection(Server_Result* res, int*& size, bigint*** Q, int number_of_clients);
	void outsource_db(string & poly_ID);
	string update(bigint elem, string delete_insert, bigint & label, string id);
	CompPerm_Request * gen_compPerm_req(byte (& tmp_key_)[AES::DEFAULT_KEYLENGTH], byte (& tmp_iv_)[AES::BLOCKSIZE]);
	
	//ning
	bigint* get_bk_ning();
	void outsource_db_ning(string & poly_ID);
	unsigned char get_bloomning();
	int  get_something();
	int bfsize_ning;
	int get_bfsize_ning();
	bigint* coefficients,  **q_ning;
	bigint** blind_coefficients, **keys_blind_coefficients;
	bigint** get_blind_coefficients();
	bigint** get_keys_blind_coefficients();
	vector <string> find_intersection_ning(bigint** q_all,Server_Result* res, int*& size, bigint*** Q, int number_of_clients, bigint*** q_ning);
	int key_A, key_B;
	bigint** find_q_all_ning(Server_Result* res, int*& size, bigint*** Q, int number_of_clients, bigint*** q_ning);
	int count_ning;
	void set_tree_q_A_ning_no_child();
	bigint **tree_q_A_ning;
	//ning

private:
	int gen_binIndx(bigint elem, int table_size);
	int* find_matched_bins(bigint k_1, bigint k_2, int size);
	int* PR_shuffle(int*elem, int size, bigint seed);
	int* find_matches(int* a, int* b, int size);
	bloom_filter convert_bigint_to_BF(bigint a, bloom_parameters parameters);
	bigint* findroots(bigint *coeff, int coeff_size, int& number_of_roots, bigint pubmoduli, ZZ_pX P);
	bigint* blind_BFs(bigint* bf, int bf_size, bigint BF_key, int bit_size, bigint pr_moduli);
	bigint* unblind_BF_(bigint BF, int indx, byte* BF_key, byte* BF_iv, bigint pr_moduli);
	bigint* blind_BF_(bigint bf, int indx, byte* BF_key, byte* BF_iv, bigint pr_moduli);
	bigint* blind_BFs_(bigint* bf, int bf_size, byte* BF_key, byte* BF_iv, bigint pr_moduli);
	bigint* unblind_BFs_(bigint* BF,int bf_size, byte* BF_key, byte* BF_iv, bigint pr_moduli);
  	bigint* gen_BF_PRN_(int indx , int counter, byte* BF_key, byte* BF_iv);
	bigint* blind_BF(bigint BF, int  indx, bigint BFkey, bigint BF_counterkey, int bit_size, bigint pr_moduli);
	bigint*	unblind_BF(bigint blinded_BF, int indx , bigint BFkey, bigint BF_ck, int bit_size, bigint pr_moduli);
  	bigint* check_vals_in_BF(bigint* vals, int val_size, bigint bf, bloom_parameters parameters, int &counter);
	bigint* assing_BFs2HT(Hashtable HT, int NoElem_in_bucket, int table_size, bloom_parameters parameters);
	bigint* gen_labels (int num, bigint seed);
	bigint* PR_shuffle(bigint* elem, int size, bigint seed);
	bigint* convert_BF_to_bigint(bloom_filter filter);
	bigint* unblind_BFs(bigint* BF, int BF_size, bigint BF_key, bigint BF_counter_key, int bit_size, bigint pr_moduli_);
	bigint** regen_bl_factors(bigint seed_, bigint ck, int* counter);
	bigint** gen_map (int size, bigint seed1, bigint seed2);
	bigint** gen_map_ (int size, byte* label_key1, byte* label_iv1, byte* label_key2, byte* label_iv2, int byte_);
	bigint** R_shuffle(bigint** elem, int size);
	bigint** PR_shuffle_bins(bigint** bins, int size, bigint seed_);
	bigint** blind_shuffled_bl(bigint**s_bl, int table_size, bigint seed_, bigint pubmoduli_);
	bigint** blind_shuffled_bl_(bigint** s_bl, int table_size, byte *key, byte* iv, int key_size, int num_of_PRNs_, int byte_, bigint pubmoduli_);
	bigint** combine_permuted_bins(bigint**& v_a, bigint**& v_b, bigint**& a, int v_size, int xpoint_size_, bigint pk_1, bigint pk_2, bigint pubmoduli);
	void get_pubModuli();
	void get_xpoints(int& size);
	void get_tablesize();
	void get_NoElem_in_bucket();
	void get_pubModuli_bitsize();
	Polynomial* PR_shuffle_poly(Polynomial* pol, int size, bigint seed);
	//Variables
	bigint  label_key, shuffle_key, BF_key, pr_moduli, *elem, pubmoduli, * xpoints;
	int key_size, elem_size, xpoint_size, *counter, table_size, pub_moduli_bitsize, labels_bit_size, pr_moduli_bitsize, NoElem_in_bucket;
	string outpoly_ID;
	Server* serv;
	bloom_parameters bf_parameters;
  	byte seed_[AES::DEFAULT_KEYLENGTH];
  	byte iv[AES::BLOCKSIZE];
  	byte BF_key_[AES::DEFAULT_KEYLENGTH];
  	byte BF_iv[AES::BLOCKSIZE];
	byte label_key_[AES::DEFAULT_KEYLENGTH];
  	byte label_iv[AES::BLOCKSIZE];
	unordered_map <string, int> xPoint_map;

	//ning
	bigint * bk_ning; //ning
	//bigint  shuffle_key_ning;
	int** assing_BFs2HT_ning(Hashtable HT, int NoElem_in_bucket, int table_size, bloom_parameters parameters);
	bigint* assing_BFs2HT_ning_bigint(Hashtable HT, int NoElem_in_bucket, int table_size, bloom_parameters parameters);
	unsigned char bloomtable_ning;//check the output
	int bf_ning[5780]; //pr_moduli_bitsize
	bigint* convert_BF_to_char_ning(bloom_filter filter);
	int** PR_shuffle_ning(int** elem, int size, bigint seed);
	int* convert_BF_to_int_ning(bloom_filter filter);
	bloom_filter shuffle_one_BF(bloom_filter filter);
	bigint* PR_shuffle_ning(bigint* elem, int size, bigint seed);
	Polynomial* PR_shuffle_poly_ning(Polynomial* pol, int size, bigint seed);
	
	//ning
};
