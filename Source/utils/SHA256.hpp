/* GRChombo - self-contained FIPS 180-4 SHA-256, no external dependency. */
#ifndef SHA256_HPP_
#define SHA256_HPP_

#include <array>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

namespace SHA256
{
inline std::uint32_t rotate(std::uint32_t x, unsigned n)
{
    return (x >> n) | (x << (32 - n));
}
inline std::string digest(std::string_view bytes)
{
    constexpr std::uint32_t k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
        0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
        0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
        0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
        0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
    std::array<std::uint32_t, 8> hash = {0x6a09e667, 0xbb67ae85, 0x3c6ef372,
                                         0xa54ff53a, 0x510e527f, 0x9b05688c,
                                         0x1f83d9ab, 0x5be0cd19};
    const std::uint64_t bits = std::uint64_t(bytes.size()) * 8;
    const std::size_t padded =
        (bytes.size() / 64 + (bytes.size() % 64 < 56 ? 1 : 2)) * 64;
    for (std::size_t offset = 0; offset < padded; offset += 64)
    {
        std::uint32_t w[64] = {};
        for (unsigned j = 0; j < 64; ++j)
        {
            const auto pos = offset + j;
            unsigned b = pos < bytes.size()
                             ? static_cast<unsigned char>(bytes[pos])
                         : pos == bytes.size() ? 0x80
                                               : 0;
            if (pos >= padded - 8)
                b = (bits >> (8 * (padded - 1 - pos))) & 0xff;
            w[j / 4] |= std::uint32_t(b) << (24 - 8 * (j % 4));
        }
        for (unsigned j = 16; j < 64; ++j)
        {
            const auto x = w[j - 15], y = w[j - 2];
            w[j] = w[j - 16] + (rotate(x, 7) ^ rotate(x, 18) ^ (x >> 3)) +
                   w[j - 7] + (rotate(y, 17) ^ rotate(y, 19) ^ (y >> 10));
        }
        auto v = hash;
        for (unsigned j = 0; j < 64; ++j)
        {
            const auto a = v[0], b = v[1], c = v[2], e = v[4], f = v[5],
                       g = v[6];
            const std::uint32_t t1 =
                v[7] + (rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25)) +
                ((e & f) ^ (~e & g)) + k[j] + w[j];
            const std::uint32_t t2 =
                (rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22)) +
                ((a & b) ^ (a & c) ^ (b & c));
            v = {t1 + t2, a, b, c, v[3] + t1, e, f, g};
        }
        for (unsigned j = 0; j < 8; ++j)
            hash[j] += v[j];
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (auto word : hash)
        out << std::setw(8) << word;
    return out.str();
}
} // namespace SHA256
#endif
