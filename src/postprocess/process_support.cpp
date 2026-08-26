#include "process_support.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <system_error>

#ifdef _WIN32
#include <process.h>
#else
#include <spawn.h>
#include <sys/wait.h>
extern char** environ;
#endif

namespace cetsp::postprocess::detail {
namespace {

std::atomic<std::uint64_t> temporaryCounter = 0;

} // namespace

TemporaryDirectory::TemporaryDirectory(std::string_view prefix) {
    const auto timestamp = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    const std::filesystem::path base = std::filesystem::temp_directory_path();

    for (std::uint64_t attempt = 0; attempt < 100; ++attempt) {
        const std::uint64_t suffix = timestamp ^
            temporaryCounter.fetch_add(1, std::memory_order_relaxed) ^ attempt;
        path_ = base / (std::string(prefix) + "-" + std::to_string(suffix));
        std::error_code error;
        if (std::filesystem::create_directory(path_, error)) {
            return;
        }
        if (error && error != std::errc::file_exists) {
            throw std::filesystem::filesystem_error(
                "could not create postprocessor temporary directory",
                path_, error);
        }
    }
    throw std::runtime_error(
        "could not create a unique postprocessor temporary directory");
}

TemporaryDirectory::~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
}

int runProcess(
    const std::filesystem::path& executable,
    const std::vector<std::string>& arguments) {
    if (!std::filesystem::is_regular_file(executable)) {
        throw std::invalid_argument(
            "postprocessor executable does not exist: " +
            executable.string());
    }

#ifdef _WIN32
    std::vector<std::wstring> wideArguments;
    wideArguments.reserve(arguments.size() + 1);
    wideArguments.push_back(executable.wstring());
    for (const std::string& argument : arguments) {
        wideArguments.push_back(
            std::filesystem::path(argument).wstring());
    }

    std::vector<const wchar_t*> argv;
    argv.reserve(wideArguments.size() + 1);
    for (const std::wstring& argument : wideArguments) {
        argv.push_back(argument.c_str());
    }
    argv.push_back(nullptr);
    return static_cast<int>(_wspawnv(
        _P_WAIT, executable.c_str(), argv.data()));
#else
    std::vector<std::string> ownedArguments;
    ownedArguments.reserve(arguments.size() + 1);
    ownedArguments.push_back(executable.string());
    ownedArguments.insert(
        ownedArguments.end(), arguments.begin(), arguments.end());

    std::vector<char*> argv;
    argv.reserve(ownedArguments.size() + 1);
    for (std::string& argument : ownedArguments) {
        argv.push_back(argument.data());
    }
    argv.push_back(nullptr);

    pid_t process = 0;
    const int spawnResult = posix_spawn(
        &process, executable.c_str(), nullptr, nullptr,
        argv.data(), environ);
    if (spawnResult != 0) {
        return spawnResult;
    }

    int status = 0;
    if (waitpid(process, &status, 0) < 0) {
        return -1;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#endif
}

} // namespace cetsp::postprocess::detail
