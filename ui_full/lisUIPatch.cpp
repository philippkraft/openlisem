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
    QNetworkAccessManager manager;
    QEventLoop loop;
    QNetworkReply *reply = manager.get(QNetworkRequest(QUrl("https://raw.githubusercontent.com/vjetten/openlisem/main_C/include/version.h")));
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();


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
        if (latestPart > currentPart) {
            return true;
        } else if (latestPart < currentPart) {
            return false;
        }
    }
    return latestParts.size() > currentParts.size();
    //return currentVersion != latestVersion;
}
//-------------------------------------------------------------------------------------
QString lisemqt::getLatestVersionFromGitHub()
{
    QNetworkAccessManager manager;
    QEventLoop loop;
    QNetworkReply *reply = manager.get(QNetworkRequest(QUrl("https://raw.githubusercontent.com/vjetten/openlisem/main_C/include/version.h")));
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
    QMessageBox msg;
    QString currentVersion = VERSIONNR;
    QString latestVersion = getLatestVersionFromGitHub();

    if (!latestVersion.isEmpty() && isNewVersionAvailable(currentVersion, latestVersion)) {

        msg.setText("A new version is available (" + latestVersion + "). Please download the latest version.");

        downloadPatch();

    } else {
        if (latestVersion.isEmpty()) {
            // Handle offline scenario
            //qDebug() << "Cannot check updates online.";
        } else {
            //msg.setText("Up to Date: \nYou are using the latest version (" + currentVersion + ").");
        }
    }

    // Create a timer to close the message box after 3 seconds
    QTimer::singleShot(3000, &msg, &QMessageBox::accept);
    msg.exec();
}



