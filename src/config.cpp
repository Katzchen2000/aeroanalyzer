#include "aeroanalyzer/config.h"
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

namespace aero {

static std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

bool Config::load(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        std::size_t h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        std::size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = trim(line.substr(0, eq));
        std::string v = trim(line.substr(eq + 1));
        if (!k.empty()) kv_[k] = v;
    }
    return true;
}

double Config::getd(const std::string& key, double def) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) return def;
    try { return std::stod(it->second); } catch (...) { return def; }
}

int Config::geti(const std::string& key, int def) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) return def;
    try { return std::stoi(it->second); } catch (...) { return def; }
}

bool Config::getb(const std::string& key, bool def) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) return def;
    std::string v = it->second;
    std::transform(v.begin(), v.end(), v.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return v == "1" || v == "true" || v == "yes" || v == "on";
}

std::string Config::gets(const std::string& key, const std::string& def) const {
    auto it = kv_.find(key);
    return it == kv_.end() ? def : it->second;
}

}  // namespace aero
