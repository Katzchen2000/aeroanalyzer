// config.h — minimal dependency-free key=value config (no YAML dep).
//
// Format: "key = value" per line, '#' starts a comment. All the magic numbers
// you'll end up tuning (masses, SM band, hinge limit, GA sizes) live in
// config/baseline.cfg so you don't recompile to retune.
#pragma once
#include <string>
#include <map>

namespace aero {

class Config {
public:
    // Loads a file if it exists; missing keys fall back to defaults at the
    // call site. Returns false if the path was given but could not be opened.
    bool load(const std::string& path);

    double getd(const std::string& key, double def) const;
    int    geti(const std::string& key, int def) const;
    bool   getb(const std::string& key, bool def) const;
    std::string gets(const std::string& key, const std::string& def) const;

    // Programmatic override (tests, CLI overrides).
    void set(const std::string& key, const std::string& val) { kv_[key] = val; }

private:
    std::map<std::string, std::string> kv_;
};

}  // namespace aero
