//
// Created by Jonathan Ebrahimian on 5/14/21.
//

#ifndef EPPSWARM_CREDENTIALS_H
#define EPPSWARM_CREDENTIALS_H
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <fstream>

class credentials {
    private:
        Aws::Auth::AWSCredentials aws_credentials;
        static credentials *instance;
        credentials() {
            std::ifstream inFile;
            std::cout << "Please provide path to AWS credentials" << std::endl;
            std::string path;
            getline (std::cin, path);
            if(path == ""){
                path = "../.env";
            }
            inFile.open(path);
            if(!inFile.is_open()){
                std::cout << "Error opening credentials file." << std::endl;
                exit(1);
            }
            std::string aws_access;
            std::string aws_secret;
            //! add error handling
            getline (inFile, aws_access);
            getline (inFile, aws_secret);
            aws_credentials.SetAWSAccessKeyId(aws_access.c_str());
            aws_credentials.SetAWSSecretKey(aws_secret.c_str());
        }

    public:
        static credentials *getInstance () {
            if(!instance){
                instance = new credentials;
            }
            return instance;
        }
        Aws::Auth::AWSCredentials getData(){
            return this->aws_credentials;
        }

};

credentials *credentials::instance = 0;
#endif //EPPSWARM_CREDENTIALS_H
