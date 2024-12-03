/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include "lisemqt.h"

#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QSslSocket>


// NOTE: on windows for a working version, openssl must be installed in msys. Then the file
// "qopensslbackend.dll" must be installed in the subdirectory tls where the exe is

//-------------------------------------------------------------------------------------

void lisemqt::downloadPatch()
{
    QString patchName = QString("openLisemSetup_%1.exe").arg(VERSIONNR);
    QString urltxt = QString("https://github.com/vjetten/openlisem/releases/download/lisem_bin/%1?raw=true").arg(patchName);

    // Ask the user if they want to download and install the patch
    QMessageBox::StandardButton replyButton;
    replyButton = QMessageBox::question(nullptr, "New LISEM version",
                                        QString("A newer LISEM version is available (%1).\nDo you want to download and install the new version?").arg(VERSIONNR),
                                        QMessageBox::Yes | QMessageBox::No);

    if (replyButton == QMessageBox::Yes) {
        QString downloadsDir = QStandardPaths::writableLocation(QStandardPaths::DownloadLocation);
        QString filePath = QDir(downloadsDir).filePath(patchName);

        // Create a file to save the patch
        QFile *file = new QFile(filePath);
        if (!file->open(QIODevice::WriteOnly)) {
            qDebug() << "Error: Unable to open file for writing.";
            delete file;
            return;
        }

        // Show download progress
        QProgressDialog *progressDialog = new QProgressDialog("Downloading patch...", "Cancel", 0, 100, nullptr);
        progressDialog->setWindowModality(Qt::WindowModal);
        progressDialog->setMinimumDuration(0);

        // Proceed with the download
        QNetworkReply *reply = manager->get(QNetworkRequest(QUrl(urltxt)));

        QObject::connect(reply, &QNetworkReply::downloadProgress, [=](qint64 bytesReceived, qint64 bytesTotal) mutable {
            if (bytesTotal > 0) {
                progressDialog->setValue(static_cast<int>((bytesReceived * 100) / bytesTotal));
                if (bytesReceived == bytesTotal) {
                    progressDialog->setLabelText("Finalizing...");
                }
            }
        });

        QObject::connect(progressDialog, &QProgressDialog::canceled, [=]() mutable {
            reply->abort();
            file->close();
            file->remove();
            QMessageBox::information(nullptr, "Download Cancelled", "The download has been cancelled.");
            progressDialog->deleteLater();
        });

        QObject::connect(reply, &QNetworkReply::readyRead, [=]() mutable {
            file->write(reply->readAll());
        });

        QObject::connect(reply, &QNetworkReply::finished, [=]() mutable {
            int err = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
            qDebug() << err;
            // Handle error scenario
            if (err != 200) {
                file->close();
                file->remove();
                QMessageBox::critical(nullptr, "Error downloading", QString("Cannot download file, return code GitHub server = %1").arg(err));
                progressDialog->deleteLater();
                return;
            }

            // Close the new file
            file->close();
            reply->deleteLater();
            file->deleteLater();

            if (reply->error() == QNetworkReply::NoError) {
                progressDialog->setValue(100);
                progressDialog->deleteLater();

                QMessageBox::StandardButton installButton;
                installButton = QMessageBox::question(nullptr, "Install new version",
                                                      "The new version has been downloaded successfully in your Download folder. Do you want to close Lisem and install the new version?",
                                                      QMessageBox::Yes | QMessageBox::No);

                if (installButton == QMessageBox::Yes) {
                    // Execute the downloaded file
                    QProcess::startDetached(filePath);

                    // Close the current application
                    QApplication::quit();
                }
            } else {
                // Handle download failure
                QMessageBox::critical(nullptr, "Download Failed", "Failed to download the new version.");
                qDebug() << "Error:" << reply->errorString();
                progressDialog->deleteLater();
            }
        });

        // Show the progress dialog
        progressDialog->show();
    }
}

//-------------------------------------------------------------------------------------
bool lisemqt::isNewVersionAvailable(QString &currentVersion, QString &latestVersion)
{
    // Assuming version strings are in the format "major.minor.patch"
    QStringList currentParts = currentVersion.split(".");
    QStringList latestParts = latestVersion.split(".");

    for (int i = 0; i < qMin(currentParts.size(), latestParts.size()); ++i) {
        int currentPart = currentParts.at(i).toInt();
        int latestPart = latestParts.at(i).toInt();

        if (currentPart > latestPart) {
            return true;
        } else if (latestPart < currentPart) {
            return false;
        }
    }

    return currentParts.size() > latestParts.size(); // does this work as default 4.8 versus 4.7.3
}
//-------------------------------------------------------------------------------------
QString lisemqt::getLatestVersionFromGitHub()
{
    QEventLoop loop;
    manager = new QNetworkAccessManager();
    QNetworkReply *reply = manager->get(QNetworkRequest(QUrl("https://raw.githubusercontent.com/vjetten/openlisem/main_C/include/version.h")));
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();

    QString latestVersion;
    if (reply->error() == QNetworkReply::NoError) {
        QByteArray response = reply->readAll();
        QString content(response);
        QRegularExpression re(R"#(#define VERSIONNR "([^"]+)")#");
        QRegularExpressionMatch match = re.match(content);
        if (match.hasMatch()) {
            latestVersion = match.captured(1);
        }
    }  else {
        // Handle the network error silently
        qDebug() << "Network error: " << reply->errorString();
    }
    reply->deleteLater();

    return latestVersion;
}
//-------------------------------------------------------------------------------------

void lisemqt::CheckVersion()
{
    QString currentVersion = VERSIONNR;
    QString latestVersion = getLatestVersionFromGitHub();

    qDebug() << currentVersion << latestVersion;

    if (!latestVersion.isEmpty() && isNewVersionAvailable(currentVersion, latestVersion)) {

        downloadPatch();

    } else {
        if (latestVersion.isEmpty()) {
            // Handle offline scenario
            //qDebug() << "Cannot check updates online.";
        } else {
            //msg.setText("Up to Date: \nYou are using the latest version (" + currentVersion + ").");
            //QTimer::singleShot(3000, &msg, &QMessageBox::accept);
            //msg.exec();
        }
    }
}



