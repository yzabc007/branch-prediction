#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

#include "utils.h"
#include "tracer.h"
#include <bitset>
#include <cmath>
#include <cstdlib>

// 3.088

#define LOG_BASE 13 //13
#define LOG_GLOBAL (LOG_BASE - 1) //1
#define N_BANKS 4 //4
#define CTR_BITS 3 //3
#define TAG_BITS 11 // 11

#define MAX_LENGTH 131 // 131
#define MIN_LENGTH 3   // 3

typedef bitset<MAX_LENGTH> history_type;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

struct base_entry {
    int pred;
    base_entry(): pred(0) {}
};

struct global_entry {
    int ctr, tag, ubit;
    global_entry(): tag(0), ubit(random() & 3) {
		ctr = (random() & ((1 << CTR_BITS) - 1)) - (1 << (CTR_BITS - 1));
	}
};

struct folded_history {
    unsigned hash;
	int MOD, ORIGIN_LEN, COMPRESSED_LEN;

    void create(int origin_len, int compressed_len) {
        hash = 0;
        ORIGIN_LEN = origin_len;
        COMPRESSED_LEN = compressed_len;
        MOD = ORIGIN_LEN % COMPRESSED_LEN;
    }

    void update(history_type h) {
        hash = (hash << 1) | h[0];
        hash ^= h[ORIGIN_LEN] << MOD;
        hash ^= (hash >> COMPRESSED_LEN);
        hash &= (1 << COMPRESSED_LEN) - 1;
    }
};

class PREDICTOR{

  // The state is defined for Gshare, change for your design

private:
    global_entry global_table[N_BANKS][1 << LOG_GLOBAL];
    base_entry base_table[1 << LOG_BASE];

    folded_history comp_hist_i[N_BANKS], comp_hist_t[2][N_BANKS];
    history_type global_history;

    int path_history;
    int G_INDEX[N_BANKS], B_INDEX;
    int lens[N_BANKS];
    int bank, alt_bank;

    int pred_store;
    bool alt_pred;


    int get_b_index(UINT32 PC) {
        return PC & ((1 << LOG_BASE) - 1);
    }

    int get_g_index(UINT32 PC, int bank) {
        int index = PC ^
                    (PC >> ((LOG_GLOBAL - N_BANKS + bank + 1))) ^
                    comp_hist_i[bank].hash;
        if (lens[bank] >= 16)
            index ^= mix_func(path_history, 16, bank);
        else
            index ^= mix_func(path_history, lens[bank], bank);
        return index & ((1 << LOG_GLOBAL) - 1);
    }

    int g_tag(UINT32 PC, int bank) {
        int temp_tag = PC ^ comp_hist_t[0][bank].hash ^ (comp_hist_t[1][bank].hash << 1);
        return temp_tag & ((1 << (TAG_BITS - ((bank + (N_BANKS & 1)) / 2))) - 1);
    }

    int mix_func(int hist, int size, int bank) {
        hist = hist & ((1 << size) - 1);
        int temp_2 = hist >> LOG_GLOBAL;
        temp_2 = ((temp_2 << bank) & ((1 << LOG_GLOBAL) - 1)) + (temp_2 >> (LOG_GLOBAL - bank));
        int temp_1 = hist & ((1 << LOG_GLOBAL) - 1);
        hist = temp_1 ^ temp_2;
        return ((hist << bank) & ((1 << LOG_GLOBAL) - 1)) + (hist >> (LOG_GLOBAL - bank));
    }

    void update_base(UINT32 PC, bool taken) {
        if (taken) {
            if (base_table[B_INDEX].pred < 1)
                base_table[B_INDEX].pred ++;
        }
        else {
            if (base_table[B_INDEX].pred > - 2)
                base_table[B_INDEX].pred --;
        }
    }

    void update_ctr(int &ctr, bool taken, int bits) {
        if (taken) {
            if (ctr < ((1 << (bits - 1)) - 1))
                ctr ++;
        }
        else {
            if (ctr > - (1 << (bits - 1)))
                ctr --;
        }
    }

    void alloc_new_hist(bool taken, UINT32 PC) {
        int minu = 3, index = 0;
        for (int i = 0; i < bank; i ++) {
            if (global_table[i][G_INDEX[i]].ubit < minu) {
                minu = global_table[i][G_INDEX[i]].ubit;
                index = i;
            }
        }
        if (minu > 0) {
            for (int i = 0; i < bank; i ++) {
                global_table[i][G_INDEX[i]].ubit --;
            }
        }
        else {
            global_table[index][G_INDEX[index]].ctr = taken ? 0 : -1;
            global_table[index][G_INDEX[index]].tag = g_tag(PC, index);
            global_table[index][G_INDEX[index]].ubit = 0;
        }
    }

    void update_hist(bool taken, UINT32 PC) {
        path_history = (path_history << 1) + (PC & 1);
        path_history &= (1 << 10) - 1;
        global_history <<= 1;
        if (taken)
            global_history |= (history_type) 1;
        for (int i = 0; i < N_BANKS; i ++) {
            comp_hist_t[0][i].update(global_history);
            comp_hist_t[1][i].update(global_history);
            comp_hist_i[i].update(global_history);
        }
    }

    void calc_index(UINT32 PC) {
        B_INDEX = get_b_index(PC);
        for (int i = 0; i < N_BANKS; i ++) {
            G_INDEX[i] = get_g_index(PC, i);
        }
    }

    void linear_search(UINT32 PC) {
        bank = alt_bank = N_BANKS;
        for (int i = 0; i < N_BANKS; i ++) {
            if (global_table[i][G_INDEX[i]].tag == g_tag(PC, i)) {
                bank = i;
                break;
            }
        }
        for (int i = bank + 1; i < N_BANKS; i ++) {
            if (global_table[i][G_INDEX[i]].tag == g_tag(PC, i)) {
                alt_bank = i;
                break;
            }
        }
    }

    bool get_base_pred(UINT32 PC) {
        return base_table[B_INDEX].pred >= 0;
    }


public:

  // The interface to the four functions below CAN NOT be changed

    PREDICTOR(void) {
        lens[0] = MAX_LENGTH - 1;
        lens[N_BANKS - 1] = MIN_LENGTH;
        for (int i = 1; i < N_BANKS - 1; i ++) {
            double temp = pow((double)(MAX_LENGTH - 1) / MIN_LENGTH, (double)i / (N_BANKS - 1));
            lens[N_BANKS - i - 1] = (int) (MIN_LENGTH * temp + 0.5);
        }
		int storage = 0;
        for (int i = 0; i < N_BANKS; i ++) {
            comp_hist_i[i].create(lens[i], LOG_GLOBAL);
            comp_hist_t[0][i].create(comp_hist_i[i].ORIGIN_LEN, TAG_BITS - ((i + (N_BANKS & 1)) / 2));
            comp_hist_t[1][i].create(comp_hist_i[i].ORIGIN_LEN, TAG_BITS - ((i + (N_BANKS & 1)) / 2) - 1);
			storage += (1 << LOG_GLOBAL) * (5 + TAG_BITS - ((i + (N_BANKS & 1)) / 2));
        }
		storage += (1 << LOG_BASE) + (1 << (LOG_BASE - 2));
		printf("storage size = %dKB\n", storage / 8 / 1024);
    }

    bool GetPrediction(UINT32 PC) {
        calc_index(PC);
        linear_search(PC);
        if (bank == N_BANKS)
            return alt_pred = get_base_pred(PC);
        if (alt_bank == N_BANKS)
            alt_pred = get_base_pred(PC);
        else
            alt_pred = (global_table[alt_bank][G_INDEX[alt_bank]].ctr >= 0);
        if (pred_store < 0 ||
                abs(2 * global_table[bank][G_INDEX[bank]].ctr + 1) != 1 ||
                global_table[bank][G_INDEX[bank]].ubit != 0)
            return global_table[bank][G_INDEX[bank]].ctr >= 0;
        return alt_pred;
    }

    void UpdatePredictor(UINT32 PC, bool resolveDir, bool predDir, UINT32 branchTarget) {
        bool alloc = (predDir != resolveDir) & (bank > 0);
        if (bank < N_BANKS) {
            bool loc_taken = global_table[bank][G_INDEX[bank]].ctr >= 0;
            bool pseudo_new_alloc = (abs(2 * global_table[bank][G_INDEX[bank]].ctr + 1) == 1) &&
                                    (global_table[bank][G_INDEX[bank]].ubit == 0);
            if (pseudo_new_alloc) {
                if (loc_taken == resolveDir)
                    alloc = false;
                if (loc_taken != alt_pred) {
                    if (alt_pred == resolveDir) {
                        if (pred_store < 7)
                            pred_store ++;
                    }
                    else {
                        if (pred_store > -8)
                            pred_store --;
                    }
                }
            }
        }
        if (alloc)
            alloc_new_hist(resolveDir, PC);
        if (bank == N_BANKS)
            update_base(PC, resolveDir);
        else
            update_ctr(global_table[bank][G_INDEX[bank]].ctr, resolveDir, CTR_BITS);
        if (predDir != alt_pred && bank < N_BANKS) {
            if (predDir == resolveDir) {
                if (global_table[bank][G_INDEX[bank]].ubit < 3)
                    global_table[bank][G_INDEX[bank]].ubit ++;
            }
            else {
                if (global_table[bank][G_INDEX[bank]].ubit > 0)
                    global_table[bank][G_INDEX[bank]].ubit --;
            }
        }
        update_hist(resolveDir, PC);
    }

    void TrackOtherInst(UINT32 PC, OpType opType, UINT32 branchTarget) {}

  // Contestants can define their own functions below

};



/***********************************************************/
#endif

