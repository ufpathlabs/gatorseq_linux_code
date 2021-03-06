CREATE TABLE gatorseq_processing
(
SAMPLE_DIR_PATH varchar(255) NOT NULL,
STATUS varchar(255) NOT NULL,
TIME_STAMP varchar(255),

MESSAGE varchar(255),
PLMO_Number varchar(255),

Test_Product_Profile varchar(255),
Test_Product_Code varchar(255),
Diagnosis varchar(255),

Primary_Tumor_Site varchar(255),
Pre_Filter varchar(255),
Report_Template varchar(255),

QCIType varchar(255),
Treatments_Policy varchar(255),
Reporting_Method varchar(255),

QCI_Upload_Message varchar(255),
QCI_Download_Message varchar(255),
EPIC_Upload_Message varchar(255),

QCI_Re_Run varchar(255),
EPIC_Re_Run varchar(255),
Perc_Target_Cells varchar(255),
Perc_Tumor varchar(255)

);

CREATE TABLE gatorseq_run
(
SAMPLE_DIR_PATH varchar(255) NOT NULL,
AccessionId varchar(255)

);

CREATE TABLE illumina_basespace_run
(
SAMPLE_NAME varchar(255) NOT NULL,
PROJECT_NAME varchar(255) NOT NULL,
GENDER varchar(255) DEFAULT 'UNKNOWN',

APP_ID varchar(255),
APP_NAME varchar(255),
APP_VERSION_SLUG varchar(255),

APP_SESSION_ID varchar(255),
APP_SESSION_STATUS varchar(255),
DATE_CREATED varchar(255),
DATE_COMPLETED varchar(255),
TOTAL_SIZE varchar(255),
RUNNING_DURATION varchar(255)

);

CREATE TABLE trusight_sample_run
(
SAMPLE_NAME varchar(255) NOT NULL,

GENDER varchar(255) DEFAULT 'UNKNOWN',
STATUS varchar(255)

);

