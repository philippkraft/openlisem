#include "lisemqt.h"
// #include <QCoreApplication>
// #include <QNetworkAccessManager>
// #include <QNetworkReply>
// #include <QNetworkRequest>
// #include <QUrl>
// #include <QTextStream>
// #include <QVersionNumber>
// #include <QDebug>
// #include <QTimer>
 #include <QSslSocket>

//const QString currentVersion = "1.0.0"; // Your current version



FileDownloader::FileDownloader(QUrl imageUrl, QObject *parent) :
 QObject(parent)
{
 connect(
  &m_WebCtrl, SIGNAL (finished(QNetworkReply*)),
  this, SLOT (fileDownloaded(QNetworkReply*))
 );

 QNetworkRequest request(imageUrl);
 m_WebCtrl.get(request);
}

FileDownloader::~FileDownloader() { }

void FileDownloader::fileDownloaded(QNetworkReply* pReply) {
 m_DownloadedData = pReply->readAll();
 qDebug() <<m_DownloadedData;
 //emit a signal
 pReply->deleteLater();
 emit downloaded();
}

QByteArray FileDownloader::downloadedData() const {
 return m_DownloadedData;
}

void lisemqt::loadImage()
{
    QPixmap buttonImage;
    buttonImage.loadFromData(m_pImgCtrl->downloadedData());
}


void lisemqt::handleNetworkReply() {
}

void lisemqt::checkForUpdates()
{
 /*
  *
  *   qDebug() << "hier" << QSslSocket::sslLibraryBuildVersionString();

    QUrl imageUrl("https://github.com/vjetten/openlisem/blob/main_C/version.txt");
    m_pImgCtrl = new FileDownloader(imageUrl, this);
    connect(m_pImgCtrl, SIGNAL (downloaded()), this, SLOT (loadImage()));

    if (QSslSocket::supportsSsl()) {
        qDebug() << "SSL is supported";
    } else {
        qDebug() << "SSL is not supported";
    }
   // return;

    QTimer timer;
    QNetworkAccessManager* manager = new QNetworkAccessManager();

    QNetworkRequest req;
    timer.start(30000);
    timer.setSingleShot(true);
    req.setUrl(QUrl("https://tweakers.net/fotoalbum/image/zUMHOtNT1EhROJTH63r39xqf.jpg"));
                    //"https://github.com/vjetten/openlisem/blob/main_C/version.txt"));
    QNetworkReply* rep = manager->post(req, "");
    QObject::connect(&timer, &QTimer::timeout, this, [&rep]() { rep->abort(); });

    QObject::connect(manager, &QNetworkAccessManager::finished, [](QNetworkReply* reply) {
        if (reply->error())
        {
            qDebug() << "ERROR!";
            qDebug() << reply->errorString();
        }
        else
        {
            qDebug() << reply->header(QNetworkRequest::ContentTypeHeader).toString();
            qDebug() << reply->header(QNetworkRequest::LastModifiedHeader).toDateTime().toString();
            qDebug() << reply->header(QNetworkRequest::ContentLengthHeader).toULongLong();
            qDebug() << reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
            qDebug() << reply->attribute(QNetworkRequest::HttpReasonPhraseAttribute).toString();
            QFile f("moderncpp.html");

            auto data = reply->readAll();
            qDebug() << data;
            if (f.open(QIODevice::WriteOnly))
            {
                f.write(data);
            }
            else
            {
                qDebug() << "error can not open file";
                qDebug() << reply->readAll();
            }
        }
    });

    QObject::connect(manager, &QNetworkAccessManager::finished, manager,
                     &QNetworkAccessManager::deleteLater);
    QObject::connect(manager, &QNetworkAccessManager::finished, rep, &QNetworkReply::deleteLater);
    */
}


