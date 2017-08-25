import FWCore.ParameterSet.Config as cms

process = cms.Process("SiStripGainValidator")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring('/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/005463E2-5185-E711-B62D-02163E011E52.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/006308DA-4B85-E711-A552-02163E01A7A4.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0085D4AB-5485-E711-A6B7-02163E01A5FF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/00892628-5085-E711-AF5B-02163E019D12.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/00CEF1FE-5485-E711-AF58-02163E0134B5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/020C26EF-5485-E711-ADF1-02163E0118FB.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/02A6D64D-4785-E711-910B-02163E01A4CA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/045CEA8B-4785-E711-B417-02163E01A354.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/047D3F2F-4285-E711-B5AB-02163E01A609.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/04843443-5085-E711-AC7A-02163E01A2C3.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/048E02C1-5485-E711-BB2C-02163E019DB7.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/04908829-5085-E711-A5FB-02163E012510.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/06475D27-4285-E711-A028-02163E01A33D.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/068E7BE9-5185-E711-86A8-02163E0143FD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/06B44651-5885-E711-91F4-02163E01A4E0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/06D6BD46-4285-E711-908B-02163E01A304.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/06E8B8D9-5185-E711-9FD2-02163E01A4E0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/080BE2CC-5185-E711-98B3-02163E019DD0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0A4B3EB3-5485-E711-81AA-02163E0144E1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0C2D6D55-4A85-E711-85B3-02163E0144DD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0C372DBB-5185-E711-BE79-02163E01A1E9.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0C791B54-4785-E711-8F2C-02163E01200E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0C8C592B-4285-E711-A36B-02163E01A37A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0CAB8FF5-4F85-E711-922C-02163E01A5E7.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0E044D31-4285-E711-914F-02163E019DD0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0E5EC477-5285-E711-AD82-02163E01A2C3.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/0EB41D42-4285-E711-A2B8-02163E0135C8.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1017B3CA-4B85-E711-B9BE-02163E013483.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/107D0927-5685-E711-9BA9-02163E01A6C5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/10902303-4385-E711-B36F-02163E01A1D3.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/10A07841-4A85-E711-AD96-02163E019B45.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/10ABA6D8-5185-E711-BD84-02163E01A4C2.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/128F3230-4285-E711-A9F9-02163E01A479.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/12A5140F-4385-E711-84A4-02163E0141EA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1459A852-4285-E711-A31F-02163E01201B.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/14F3883E-4A85-E711-B1E6-02163E01A2C0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/167B484F-4A85-E711-9EFA-02163E01A43A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/16D64142-5885-E711-BA08-02163E01A523.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/180E117C-4785-E711-BF2F-02163E011AFE.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/18129F2B-4285-E711-A876-02163E01A4CF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/181DC31F-4285-E711-8639-02163E019C10.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/18956070-3B85-E711-B202-02163E01A609.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1898A670-4785-E711-B05F-02163E019DBF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/18B2D65C-5785-E711-90E2-02163E019D0A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/18FED578-3B85-E711-96CD-02163E01A5E1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1A585050-4285-E711-BE95-02163E011838.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1A8C8F5A-4785-E711-9D46-02163E01A5D1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1ACF7C61-5685-E711-AA64-02163E019B42.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1C2DB840-4A85-E711-8EBD-02163E0142EC.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1C4E86E5-4285-E711-99BE-02163E01A1E4.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1C5EF086-4A85-E711-BF83-02163E01381D.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1C76FE3B-4285-E711-BDE3-02163E01A42C.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1C94EF21-4385-E711-901E-02163E013630.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1E16331B-5085-E711-9DF7-02163E011DE8.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1E799650-4A85-E711-A968-02163E01A518.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1E8EE337-4285-E711-83CC-02163E01391D.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1EA768D2-5485-E711-880F-02163E0137F6.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1EB67433-4285-E711-A648-02163E0145DA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1EBE2A4A-4A85-E711-8CB2-02163E01A4F5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/1ED1C63E-4285-E711-98A7-02163E011E55.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/202CDE6A-4285-E711-B07A-02163E0135F6.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/203A1A26-4285-E711-A434-02163E01A43A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/20860C33-4285-E711-8FBB-02163E0143D6.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/209669F4-5485-E711-9AEC-02163E0144AC.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/20B5D97A-5385-E711-9162-02163E01A66C.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/20B63B59-3F85-E711-B941-02163E011A01.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/20D111FF-4D85-E711-9200-02163E0138BB.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22077E01-5285-E711-81EB-02163E0117FF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2207C327-5085-E711-B96B-02163E013705.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2238D9F3-3C85-E711-8B1F-02163E01A48E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22479F7F-3B85-E711-BC27-02163E019B9F.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22869888-3B85-E711-AEDE-02163E011EF1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22CF4505-5085-E711-A850-02163E01A4AE.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22DE514A-4A85-E711-8E9F-02163E01A6CA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/22E3DE2F-5085-E711-AFDC-02163E012A34.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2402922D-4285-E711-9DC1-02163E01A48E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/240A520F-4385-E711-A37B-02163E012AFE.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2423D41A-3D85-E711-8306-02163E011B87.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/247C2820-5885-E711-A9B3-02163E01A3D2.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/24AB1F33-4285-E711-A24B-02163E0144DD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/24C104F2-3C85-E711-B3DA-02163E01A1BD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/24DEFA3A-4585-E711-85F7-02163E0146B5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/262948FE-3C85-E711-BA9F-02163E01A4E9.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/262D41F6-5185-E711-B77E-02163E0141E8.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2691C040-4285-E711-8616-02163E019BDF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/26E79C2C-4285-E711-A588-02163E01A704.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2829DF50-4785-E711-989B-02163E01451D.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2898E88C-4285-E711-97E0-02163E01189B.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/28A451FB-4F85-E711-9656-02163E01A377.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/28B710F2-4F85-E711-9390-02163E01A2B9.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/28C10F02-5085-E711-A7ED-02163E01A2D8.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2A6AC608-3D85-E711-BEF4-02163E012A7E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2A83D6A9-5485-E711-81F6-02163E019BBE.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2AE0DD45-4A85-E711-86DE-02163E0143CF.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2C05B9E6-5185-E711-96BA-02163E019CE1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2C20C8D1-4B85-E711-B717-02163E01A61E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2CBED29D-3B85-E711-AAD2-02163E019E54.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2CDFC834-4285-E711-A78A-02163E01A73A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2CE2AA8F-4D85-E711-A957-02163E019DD0.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2CECC52F-5685-E711-9ED9-02163E01A1DD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2E36D118-5285-E711-96A2-02163E011EBA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2E5A7533-4285-E711-8700-02163E019C73.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/2EDCF67C-3B85-E711-B62C-02163E01A38E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3016A01B-4985-E711-A944-02163E01366D.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/305AD31A-5585-E711-BDB2-02163E01A379.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/306D2066-4A85-E711-A893-02163E01A1E9.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/30A4AA26-5085-E711-93E4-02163E01441B.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/30CA15D6-4B85-E711-A8CF-02163E01440E.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/30DDEB69-4785-E711-9966-02163E011EBA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/321EE953-4A85-E711-A4C5-02163E01A40A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/327384F4-4F85-E711-8298-02163E0142C5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3279F5CC-4B85-E711-9274-02163E019B52.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/32D72E1F-5885-E711-8B01-02163E0142C6.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/34136A07-4285-E711-8B90-02163E01A4E9.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/342B6626-5085-E711-8615-02163E019E58.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3470E477-3B85-E711-893C-02163E019B86.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3483356C-4785-E711-B796-02163E0142C5.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/34BABAAF-4285-E711-87E0-02163E011CDB.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/34F9ACC5-4B85-E711-871D-02163E01A23F.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/36239141-4285-E711-BE7F-02163E0133BA.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/362D5642-4A85-E711-BD60-02163E01473A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3698E0E0-5185-E711-8DDE-02163E011D31.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/36C7E842-4A85-E711-8017-02163E011A94.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/382FE732-4285-E711-B5FA-02163E01A2B1.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3835CF91-3B85-E711-B01F-02163E011E55.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3838BB55-4785-E711-BA43-02163E019C6B.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/38D16AD7-5485-E711-917E-02163E0146D7.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3A02274D-4A85-E711-86DD-02163E014217.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3A068C43-4A85-E711-8B76-02163E01A379.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3C01F5E9-5185-E711-8179-02163E019BFD.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3C65A931-4285-E711-9CDC-02163E01A36A.root',
                                                              '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBiasAAG-Express-v3/000/301/461/00000/3C914E3D-4285-E711-877E-02163E01A4DF.root',

# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/34E2149B-3485-E711-AF48-02163E01A332.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D6564C87-3B85-E711-BB6B-02163E019C3B.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6421B77D-3B85-E711-9A58-02163E01A5BD.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/B4A0BF91-3B85-E711-B01F-02163E011E55.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/88309B7A-3B85-E711-9097-02163E019DCF.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/323534F0-3C85-E711-B196-02163E01A3CE.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/DEF43A19-3D85-E711-86B7-02163E014138.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6CB08228-3D85-E711-8C76-02163E011A01.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/AAE49908-3D85-E711-BEF4-02163E012A7E.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D414AF2F-4285-E711-9A79-02163E019CCB.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/C448402F-4285-E711-BD41-02163E01432C.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/828BFA2D-4285-E711-B18D-02163E01440E.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7E85E420-4285-E711-8885-02163E01A4F5.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/B0FF2920-4285-E711-86EA-02163E019CD0.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/C481C723-4285-E711-88B9-02163E01A3D2.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D6CF070F-4385-E711-8CCD-02163E0141EA.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/14B185E5-4285-E711-99BE-02163E01A1E4.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/0A126D10-4385-E711-A9FB-02163E0142B1.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/96081603-4385-E711-85D2-02163E01A1D3.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/F888CC3D-4285-E711-A81F-02163E019DDD.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/4C04EC35-4285-E711-98EA-02163E01A674.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7C235249-4285-E711-A75F-02163E01262C.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/58BC536E-4285-E711-85AE-02163E01A20D.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/20036652-4285-E711-86BC-02163E01A2C3.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/22602AA2-4285-E711-A472-02163E011E2B.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/1CA79F25-4585-E711-942A-02163E0136F1.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6A67EE36-4585-E711-80F6-02163E01A5CE.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/DEA89555-4785-E711-8014-02163E019CD2.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7E940E56-4785-E711-AAEC-02163E012510.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/5E048D5A-4785-E711-9D46-02163E01A5D1.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/98CB945B-4785-E711-80E8-02163E019DF3.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/52981554-4785-E711-8F2C-02163E01200E.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/940AB340-4A85-E711-8EBD-02163E0142EC.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7A3A5042-4A85-E711-BD60-02163E01473A.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/428E1863-4A85-E711-A5D9-02163E01A7A4.root',
# '/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7A44C144-4A85-E711-BC19-02163E019C13.root', 
 )                            
                            )


from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
process.AAGFilter = triggerResultsFilter.clone(triggerConditions = cms.vstring("HLT_ZeroBias_FirstCollisionAfterAbortGap_*"),
                                               hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                               l1tResults = cms.InputTag( "" ),
                                               throw = cms.bool(False)
                                               )

###################################################################
# The track refitter
###################################################################
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'92X_dataRun2_Express_v4','')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("SiStripApvGain3Rcd"),
             tag = cms.string("SiStripGainFromParticles"),
             connect = cms.string("sqlite_file:/afs/cern.ch/user/m/musich/public/forStripDBcontacts/G2_initial2017/Gains_G2_299061_Sqlite.db")
             )
    )

###################################################################
# The BeamSpot
###################################################################
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

###################################################################
# The track refitter
###################################################################
from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone(src = cms.InputTag("ALCARECOSiStripCalMinBiasAAG"))

###################################################################
# The module
###################################################################
process.SiStripGainsValidator = cms.EDAnalyzer('SiStripGainsValidator',
                                               Tracks      = cms.InputTag("tracksRefit")
                                               )

###################################################################
# Output name
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("SiStripGainValidation.root")
                                   )

process.p = cms.Path(process.offlineBeamSpot*process.MeasurementTrackerEvent*process.tracksRefit*process.SiStripGainsValidator)
