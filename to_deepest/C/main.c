#include <stdint.h>
#include <stdio.h>
#include "data.h"    // 假設 data.c, data.h 中包含 f, g, P, omegas 等定義
#include "ntt.h"     // 假設 ntt.h, ntt.c 中實作 void ntt(int32_t *poly, int n, int prime_idx)
#include "ecrt.h"    // 假設 ecrt.h, ecrt.c 中有 void ecrt_rns_extension(int32_t *u, int n, int old_modnum, int new_modnum)

int num_of_used_moduli_fg[11] = {1,1,1,2,4,8,16,28,56,108,216};

/*
本程式碼展示 Falcon Keygen 計算 N^(k)(f), N^(k)(g) 每層的 NTT 結果的儲存策略。

主要特點：
1. 每層計算完成後 (取得 N^(k)(f), N^(k)(g) 的 NTT 結果) 都會將該層結果 append 到 temp 後方。
2. 資料順序為：
   NTT(f), NTT(g), NTT(N(f)), NTT(N(g)), NTT(N^2(f)), NTT(N^2(g)), ... , NTT(N^{10}(f)), NTT(N^{10}(g))
   即每層先放 f 的結果，再放 g 的結果，然後進入下一層。
3. 每層的結果以 coefficient-major 方式儲存，並在 RNS 下的 NTT domain。
4. 在初始 (layer=0) 層，我們只有 f,g 未做 NTT，所以要先對 f,g 作 NTT 得到 NTT(f), NTT(g) 並存入 temp。
5. 接著計算 N(f), N(g) 並進行 RNS 擴展(若需要)、NTT，再存 NTT(N(f)), NTT(N(g))，依此類推直到 N^{10}(f), N^{10}(g)。

ASCII ART 舉例：
假設初始時 (layer=0):
temp 佈局：
|---NTT(f)---|---NTT(g)---|
f,g 都是長度=1024, mod_count=1

下一層 (layer=1) 完成後：
|---NTT(f)---|---NTT(g)---|---NTT(N(f))---|---NTT(N(g))---|
N(f),N(g) 長度=512,可能 mod_count仍=1

layer=2後：
|---NTT(f)---|---NTT(g)---|---NTT(N(f))---|---NTT(N(g))---|---NTT(N^2(f))---|---NTT(N^2(g))---|
N^2(f),N^2(g) 長度=256，假設mod_count=1（若仍是1 mod）

假設在某層 mod_count增加了，則該層結果區段會更大。

最終layer=10：
|NTT(f)| NTT(g)| NTT(N(f))| NTT(N(g))| ... |NTT(N^{10}(f))| NTT(N^{10}(g))|

這些資料皆保存在 temp 不重寫先前層的結果。
*/

/*
輔助函式：計算某層結果所需的 int32_t 數量
length_k = 1024 / (2^k) 為第 k 層的多項式長度 (對N^(k)(f))
mod_count_k = num_of_used_moduli_fg[k]
f,g 各需要 length_k * mod_count_k int32_t
合計該層 (f或g) 的資料量為 length_k * mod_count_k
該層 f,g 合計為 2 * length_k * mod_count_k (如果要一次空間規劃可用，但我們在append時只需記住起始點)
*/

// 計算 N(f), N(g) 的函式 (NTT domain)
void compute_field_norm_once(int32_t *f_poly_ntt, int32_t *g_poly_ntt, int length, int mod_count) {
    int half_len = length / 2;

    // ASCII ART before compute:
    // f_poly_ntt: [f(0,mod0..mod_{M-1}), f(1,mod0..mod_{M-1}), ..., f(length-1, mod...)]
    // g_poly_ntt 同理
    //
    // After compute:
    // NTT(N(f)(i)) = f_ntt(2i)*f_ntt(2i+1)
    // NTT(N(g)(i)) = g_ntt(2i)*g_ntt(2i+1)
    //
    // 結果覆蓋到前 half_len 區段

    for (int i = 0; i < half_len; i++) {
        for (int m = 0; m < mod_count; m++) {
            int32_t x_f = f_poly_ntt[2*i*mod_count + m];
            int32_t y_f = f_poly_ntt[(2*i+1)*mod_count + m];
            int32_t mod = P[m]; // 簡化，不考慮 prime_idx mapping
            int64_t prod_f = (int64_t)x_f * y_f;
            f_poly_ntt[i*mod_count + m] = (int32_t)(prod_f % mod);

            int32_t x_g = g_poly_ntt[2*i*mod_count + m];
            int32_t y_g = g_poly_ntt[(2*i+1)*mod_count + m];
            int64_t prod_g = (int64_t)x_g * y_g;
            g_poly_ntt[i*mod_count + m] = (int32_t)(prod_g % mod);
        }
    }
}


// 主程式
void compute_all_field_norms(int32_t *temp) {
    // 初始化將 f,g 複製並做 NTT(f), NTT(g)
    // 最初 mod_count=1, length=1024
    // 將 NTT(f) 從 temp[0] 開始存放, 長度=1024, mod_count=1 => 1024 ints
    // NTT(g) 接在其後 => 1024 ints

    // ASCII ART 初始放入 temp（未NTT前）：
    // temp:
    // f: | f(0) | f(1) | ... | f(1023) |
    // g: | g(0) | g(1) | ... | g(1023) |
    // length=1024, mod_count=1

    for(int i = 0; i < 1024; i++) {
        temp[i] = f[i];
        temp[i + 1024] = g[i];
    }

    int32_t *f = temp;
    int32_t *g = temp + 1024;

    // 對 f,g 做 NTT
    ntt(f, 1024, 0); // NTT(f)
    ntt(g, 1024, 0); // NTT(g)

    // temp:
    // |---- NTT(f) ----|---- NTT(g) ----|


    int layer; 
    int offset = 0; // offset 用於指示在 temp 中 append 的位置
    int length = 1024;
    int current_mod_count = 1;


    // 此時 
    // temp[   0..1023] = NTT(f) (1 mod)
    // temp[1024..2047] = NTT(g) (1 mod)
    // offset = 0 為 NTT(f) 起點, NTT(g) 起點= 1024
    // append 位置將在 NTT(g) 之後
    // 計算 NTT(f)區段大小: 1024 * 1 = 1024
    // 計算 NTT(g)區段大小: 同 1024
    // 下一層結果 append 起點 = 2048
    int ntt_f_size = length * current_mod_count; // 1024
    int ntt_g_size = length * current_mod_count; // 1024
    int next_append_offset = ntt_f_size + ntt_g_size; // 2048

    // ASCII ART:
    // temp:
    // |---- NTT(f) ----|---- NTT(g) ----|
    // len(NTT(f))=1024, len(NTT(g))=1024
    // index range:
    // NTT(f): [0..1023]
    // NTT(g): [1024..2047]
    // next_append_offset = 2048

    // 現在開始計算 N^(1)(f), N^(1)(g)
    for (layer = 1; layer <= 10; layer++) {
        int required_mod_count = num_of_used_moduli_fg[layer];
        // 若需擴展 RNS:
        if (required_mod_count > current_mod_count) {
            // intt 
            f = temp; // NTT(f) 的起點
            g = temp + length * current_mod_count; // NTT(g) 的起點

            // ecrt_rns_extension

            // ntt

            current_mod_count = required_mod_count;
        }

        // length 已是前一層的長度, 現在還是上一層的 N^(layer-1)(f,g)
        // 我們需要對 N^(layer-1)(f,g) 做 NTT (若尚未做 or 需更新)
        // 事實上N^(layer-1)(f,g) 已在 NTT domain 中（依前面設計，每層結果都是NTT狀態保存）
        // 所以此處可能不需要再 NTT，只要確保上一層已計算完畢是 NTT domain 即可
        // 假設每層計算結束結果都是 NTT domain:
        // 那麼 N^(layer-1)(f,g) 已經在 temp 中最後 append 的區域

        // 找出上一層 N^(layer-1)(f) 的位置
        // 前面層層 append:
        // layer=0 時放 NTT(f),NTT(g)到 [0..2047]
        // layer=1 時 append NTT(N(f)),NTT(N(g)) 後面
        // 因為每層尺寸不同，需要動態計算offset
        // 這裡我們用一個簡化方式：
        // 就假設我們每層計算完，當下記住該層 f,g 的起始位置與大小。
        // （實務中可以有個表記錄每層的 offset與大小）

        // 此範例：假設我們在每層計算後將結果 append 並立即更新 offset 等
        // 上一層的 f,g 結果已經存放完成，我們可以從 next_append_offset 開始存當前層的結果。

        int prev_layer_length = length; // 上一層的 length
        int prev_layer_mod_count = current_mod_count; // 上一層的mod_count
        int prev_ntt_f_size = prev_layer_length * prev_layer_mod_count;
        int prev_ntt_g_size = prev_ntt_f_size; // f,g大小相同
        int prev_layer_total_size = prev_ntt_f_size + prev_ntt_g_size;

        // 上一層結果結尾即是下一層寫入起點:
        int prev_layer_start = next_append_offset - prev_layer_total_size;
        int f_poly_ntt_start = prev_layer_start;              // N^(layer-1)(f)的起點
        int g_poly_ntt_start = prev_layer_start + prev_ntt_f_size; // N^(layer-1)(g)的起點

        // 現在計算 N^(layer)(f),N^(layer)(g):
        // N^(layer)(f) = N^(layer-1)(f)(2i)*N^(layer-1)(f)(2i+1)
        // 同理 g
        compute_field_norm_once(temp + f_poly_ntt_start, temp + g_poly_ntt_start, prev_layer_length, prev_layer_mod_count);

        // 計算後 length 減半
        length = prev_layer_length / 2;

        // 計算本層結果大小
        int this_ntt_f_size = length * prev_layer_mod_count; 
        int this_ntt_g_size = this_ntt_f_size;
        int this_layer_total_size = this_ntt_f_size + this_ntt_g_size;

        // 現在 N^(layer)(f),N^(layer)(g) 已經在 temp[f_poly_ntt_start ...] (前 half_len 部分)
        // 但是這樣會覆蓋上一層的資料嗎？
        // 上一層結果是 append 在後面，不會覆蓋前面。
        // 但目前我們計算都是在原地修改上一層的 f,g 區段的前 half_len 區段。
        // 要將本層結果 append 到 temp 的後方去儲存，以保留各層結果。

        // STEP:
        // 1. 本層結果現在在 f_poly_ntt_start,g_poly_ntt_start處(前 half_len*mod_count)
        // 2. 我們需將這些結果複製到 next_append_offset 開始的地方（append）
        //    以保留本層結果到後面區域
        int32_t *f_src = temp + f_poly_ntt_start; 
        int32_t *g_src = temp + g_poly_ntt_start;

        int32_t *f_dest = temp + next_append_offset; // append本層 f 結果的位置
        int32_t *g_dest = f_dest + this_ntt_f_size;  // append本層 g 結果的位置

        // 複製本層 f,g 結果到 append 區域
        for (int i = 0; i < this_ntt_f_size; i++) {
            f_dest[i] = f_src[i];
        }
        for (int i = 0; i < this_ntt_g_size; i++) {
            g_dest[i] = g_src[i];
        }

        // ASCII ART:
        // temp結構 (以layer=1為例):
        // | NTT(f) 1024 ints | NTT(g) 1024 ints | NTT(N(f)) 512 ints | NTT(N(g)) 512 ints |
        // idx range:
        // NTT(f)   [0..1023]
        // NTT(g)   [1024..2047]
        // NTT(N(f))[2048..2048+512-1=2559]
        // NTT(N(g))[2560..2560+512-1=3071]

        // 更新 next_append_offset 為下層append起點
        next_append_offset += this_layer_total_size;

        // current_mod_count 在本例中若有擴展須更新
        // 若已擴展則 current_mod_count = required_mod_count;
        // 否則不變。

        // 下一層使用更新後的 length, current_mod_count 繼續
    }

    // 迄今為止：
    // temp中按照:
    // NTT(f), NTT(g),
    // NTT(N(f)), NTT(N(g)),
    // NTT(N^2(f)), NTT(N^2(g)), 
    // ...
    // NTT(N^{10}(f)), NTT(N^{10}(g))

    // 都以append方式存下來。
    // f,g 在計算 N^(k)(f,g)時不相互影響，但我們保留了所有層的 f,g 結果在 temp 中。
}

int main() {
    // 為了確保有足夠空間來存所有層的結果 (f,g)：
    // 最大層使用 216 個模數, 最小長度為1係數(到N^{10})，最大長度1024, 所有層記錄。
    // 粗略估計 max space: 每層f,g資料量合計不斷疊加，比 1024*216*2 大得多。
    // 此處給一個很大的暫定值，比如:
    static int32_t temp[1024*216*2*20]; // 實務依需求調整

    compute_all_field_norms(temp);

    // 最終 temp 中有所有層的 NTT(N^(k)(f)), NTT(N^(k)(g))。
    // 後續步驟可依需要索引相應區段使用。

    return 0;
}
