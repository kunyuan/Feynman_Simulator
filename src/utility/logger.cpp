#include <iostream>
#include <new>
#include <cstdlib>
#include <ctime>

#include "logger.h"

// Definition (and initialization) of static attributes
Logger* Logger::m_ = 0;

#ifdef LOGGER_MULTITHREAD
pthread_mutex_t Logger::lock_ = PTHREAD_MUTEX_INITIALIZER;
inline void Logger::lock()
{
    pthread_mutex_lock(&lock_);
}

inline void Logger::unlock()
{
    pthread_mutex_unlock(&lock_);
}
#else
void Logger::lock()
{
}
void Logger::unlock()
{
}
#endif

/**
 * \brief Constructor.
 * It is a private constructor, called only by getInstance() and only the
 * first time. It is called inside a lock, so lock inside this method
 * is not required.
 * It only initializes the initial time. All configuration is done inside the
 * configure() method.
 */
Logger::Logger()
    : configured_(false)
{

    gettimeofday(&initialTime_, NULL);
}

/**
 * \brief Method to configure the logger. Called by the DEBUG_CONF() macro.
 * To make implementation easier, the old stream is always closed.
 * Then, in case, it is open again in append mode.
 * @param Name of the file used for logging
 * @param Configuration (i.e., log on file and on screen on or off)
 * @param Verbosity threshold for file
 * @param Verbosity threshold for screen
 */
void Logger::configure(const std::string& outputFile,
                       const std::string& loggerName,
                       const loggerConf configuration,
                       const int fileVerbosityLevel,
                       const int screenVerbosityLevel)
{
    Logger::lock();
    loggerName_ = "[" + loggerName + "]";
    fileVerbosityLevel_ = fileVerbosityLevel;
    screenVerbosityLevel_ = screenVerbosityLevel;

    // Close the old stream, if needed
    if (configuration_ & file_on)
        out_.close();

    // Compute a new file name, if needed
    if (outputFile != logFile_) {
        logFile_ = outputFile;
    }

    // Open a new stream, if needed
    if (configuration & file_on)
        out_.open(logFile_.c_str(), std::ios::app);

    configuration_ = configuration;
    configured_ = true;

    Logger::unlock();
}

/**
 * \brief Destructor.
 * It only closes the file, if open, and cleans memory.
 */

Logger::~Logger()
{
    Logger::lock();
    if (configuration_ & file_on)
        out_.close();
    delete m_;
    Logger::unlock();
}

/**
 * \brief Method to get a reference to the object (i.e., Singleton)
 * It is a static method.
 * @return Reference to the object.
 */
Logger& Logger::getInstance()
{
    Logger::lock();
    if (m_ == 0)
        m_ = new Logger;
    Logger::unlock();
    return *m_;
}

/**
 * \brief Method used to print messages.
 * Called by the DEBUG() macro.
 * @param Priority of the message
 * @param Source file where the method has been called (set equal to __FILE__
 * 	      by the DEBUG macro)
 * @param Source line where the method has been called (set equal to __LINE__
          by the macro)
 * @param Message
 */
void Logger::print(const unsigned int verbosityLevel,
                   const std::string& path,
                   const int line,
                   const std::string& message)
{
    if (!configured_) {
        std::cout << "ERROR: Logger not configured!" << std::endl;
        return;
    }

    size_t sep = path.find_last_of("\\/");
    std::string file;
    if (sep != std::string::npos)
        file = path.substr(sep + 1, path.size() - sep - 1);

    struct timeval currentTime;
    gettimeofday(&currentTime, NULL);
    time_t _time = time(NULL);
    std::string str(ctime(&_time));

    std::ostringstream oss;
    oss << loggerName_ << LOGSTR[verbosityLevel];
    time_t currTime;
    time(&currTime);
    struct tm* currTm = localtime(&currTime);
    oss << "[" << (currTm->tm_year - 100) << "/" << currTm->tm_mon << "/" << currTm->tm_mday << " " << currTm->tm_hour << ":" << currTm->tm_min << ":" << currTm->tm_sec << "]";
    oss << "@[" << file << ":" << line << "]\n" << message << std::endl;
    std::string msg = oss.str();
    Logger::lock();
    if ((configuration_ & file_on) && (verbosityLevel >= fileVerbosityLevel_))
        out_ << msg << std::endl;
    if ((configuration_ & screen_on) && (verbosityLevel >= screenVerbosityLevel_))
        std::cout << msg << std::endl;
    Logger::unlock();
}
