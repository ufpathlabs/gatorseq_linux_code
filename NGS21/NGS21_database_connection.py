import yaml
import os
import mysql.connector
def getSQLConnection(config_path, token_path, CODE_ENV):
    config_dict=dict()
    with open(config_path, 'r') as stream:
        try:
            config_dict=yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()

    USER_NAME=os.environ['USER']

    def replace_env(strname):
        strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
        return strname

    config_token_dict=dict()
    with open(token_path, 'r') as stream:
        try:
            config_token_dict=yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit()
    
    if CODE_ENV == "DevEnv":
        MYSQL_HOST = config_dict['DEV_MYSQL_HOST']
        MYSQL_USERNAME = config_dict['DEV_MYSQL_USERNAME']
        MYSQL_DATABASE = config_dict['DEV_MYSQL_DATABASE']
        MYSQL_PASSWORD = config_token_dict['DEV_MYSQL_PASSWORD']

    if CODE_ENV == "ProdEnv":
        MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
        MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
        MYSQL_PASSWORD = config_token_dict['PROD_MYSQL_PASSWORD']
        MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']

    conn = None
    try:
        conn = mysql.connector.connect(
            host=MYSQL_HOST,
            user=MYSQL_USERNAME,
            passwd=MYSQL_PASSWORD,
            database=MYSQL_DATABASE
        )
    except:
        print(traceback.format_exc())
        print("----- unable to get connection-----")
    return conn
