---
description: Manage your CyteType reports, API keys and account in the dashboard, in your personal workspace or your organization's workspace.
---

(dashboard)=
# Dashboard

The dashboard at `https://cytetype.nygen.io/dashboard` is where you manage your reports, your API keys and your account. It needs a signed-in browser session: if you are not signed in, it sends you to the sign-in page first and brings you back afterwards. See [Sign in and API keys](../getting-started/sign-in.md).

## Open the dashboard

Open the URL in a browser, or open it from your terminal or R session:

::::{tab-set}
:::{tab-item} Python
:sync: python
```bash
cytetype dashboard
```
:::
:::{tab-item} R
:sync: r
```r
OpenCyteTypeDashboard()
```
:::
::::

Both open the dashboard of the server you signed in to, or `https://cytetype.nygen.io/dashboard` when there is no saved sign-in. They also print the URL, so you can copy it when no browser opens.

The sidebar links to **Chat**, **Reports**, **API keys** and the documentation. At the bottom are the workspace switcher (only for members of an organization with an active license) and your account.

## Workspaces

Every account has a personal workspace. If your account belongs to an organization with an active license, the sidebar also shows a workspace switcher with your organization (and your role in it) and **Personal**. The dashboard offers one organization per account: the one it joined first. The Reports, API keys and Account pages always show the selected workspace only.

| | Personal workspace | Organization workspace |
|---|---|---|
| Who has it | Every account | Members of an organization with an active license |
| Reports | Jobs submitted with your personal API keys | Jobs submitted with your organization's API keys |
| API keys | Personal keys (`cyt_p_...`) that you create | Organization keys (`cyt_...`) that admins create |
| Usage counted against | Your personal rolling limits | The organization's cluster quota |
| Visibility options | Private, Public | Private, Organization, Public |

The API key that submits a job decides its workspace. A job submitted with a personal key belongs to your personal workspace and counts toward your personal limits. A job submitted with an organization key belongs to the organization workspace and counts toward the organization's quota. Browser sign-in with `cytetype setup` or `SetupCyteTypeR()` gives you a personal key, so to run jobs under your organization, use an organization key from an admin (see [Sign in and API keys](../getting-started/sign-in.md)).

When an organization's license expires, its workspace is no longer offered in the switcher.

## Reports

The **Reports** page has two tabs. The first is the same in both workspaces:

- **My Reports**: the reports you own in this workspace, whether you submitted them or they were transferred to you.

The second depends on the workspace:

- **Shared with Me** (personal workspace): personal reports that other people shared with you by email.
- **Workspace Reports** (organization workspace): other members' reports with Organization or Public visibility, and their private reports that were shared with you.

Each row in **My Reports** shows the report name, the report ID (with a copy button), the status (Processing, Completed or Failed), the visibility, who it is shared with, the number of clusters, the creation date and the actions. The other tab shows the owner instead of the sharing details.

:::{note}
A report that belongs to an organization you are not a member of does not appear in either tab, even when it was shared with you. Open it from the link in the share email.
:::

### Open a report, Chat or target analysis

Each row has up to three buttons:

- **View** opens the report in a new tab (see [The CyteType report](../results/report.md)).
- **Chat** opens [Chat](chat.md) in a new tab with this run attached. It appears for completed reports only.
- **Target analysis** opens the [target analysis](target-analysis.md) dialog. It appears for completed reports only.

### Rename a report

In **My Reports**, click the report name, type the new name and press Enter (Esc cancels). Reports without a name show "Untitled". You can also name a report when you submit it, with a `title` key in `metadata` (see [Study context and metadata](../guides/study-context.md)). Only the owner can rename a report.

### Change the visibility

Click the edit icon next to the visibility badge, choose **Private**, **Organization** or **Public**, and click **Save**. **Organization** is offered only in the organization workspace. You can change the visibility in either direction at any time. [Visibility and sharing](sharing.md) explains what each level allows.

### Share by email

1. Tick the reports to share and click **Share**, or click the **Shared With** cell of a single report.
2. Enter an email address and click **Add** (or press Enter). Repeat for more people.
3. Click **Save Changes**.

What happens next:

- Each recipient gets view access and one email that links to every report newly shared with them in this save.
- A recipient without a CyteType account is listed as **Pending**. The invite turns into access the first time they sign in with that address.
- Sharing again with someone who already has access does not send another email.
- Recipients need an institutional email address, the same rule as for signing in.
- With several reports selected, the dialog only adds people. To remove someone, including a pending invite, open the dialog for that one report, remove the address and click **Save Changes**.

Only the owner can share a report. Owners can also share from the **Share** button on the report page.

### Transfer ownership

Transfer is available in the organization workspace only.

1. Tick the reports and click **Transfer**.
2. Enter the new owner's email address and click **Transfer All**.
3. Confirm.

The new owner must be a member of the same organization. A transfer cannot be undone: the new owner gets all owner rights (rename, visibility, sharing, deletion) and you lose them. If the report is private, you also lose access to it unless the new owner shares it with you or sets it to Organization.

### Archive a report

Tick the reports and click **Archive**. Archived reports disappear from the list until you tick **Show archived**; select them there and click **Unarchive** to bring them back.

Archiving does not change who can open the report: the link keeps working for everyone who could open it before. It removes the report from your lists, and from the report search in Chat and the job list of the [MCP connector](mcp.md) for everyone. An archived report can still be attached in Chat or connected through MCP by its full report ID.

### Delete a report

Tick the reports in **My Reports**, click **Delete** and confirm. Only the owner can delete a report.

Deletion is a soft delete:

- It cannot be undone from the dashboard. Deleted reports stay in **My Reports** under **Show deleted**, without any actions.
- Everyone loses access, including you. All shares and pending invites are removed, and anyone who opens the link sees "Report Deleted: This report has been deleted by the owner."
- Requests for its results with an API key fail with the error code `REPORT_DELETED`.
- Chat and the MCP connector can no longer use it. Chat answers and outputs that were based on it become unavailable.
- Deleted jobs still count toward your usage limits or your organization's quota (see [Usage limits and quotas](limits.md)).

If you delete a report that is still processing, it disappears from your views but the job keeps running on the server until it finishes.

## API keys

API keys let the Python and R clients submit jobs and fetch results. The **API keys** page lists the keys of the selected workspace.

### Create a key

1. Select the workspace the key is for.
2. Click **Create API key**.
3. Enter a **Key name**, for example the computer or pipeline that will use it.
4. Optionally enter a **Quota limit**: the maximum number of clusters this key can annotate in total. Leave it empty for no per-key limit. Your other limits still apply.
5. Organization admins can enter a member's address in **User email (optional)** to create the key for that member.
6. Click **Create API key**, then copy the key from the **API key created** panel.

:::{warning}
The dashboard shows a key you create here only once. If you lose it, disable it and create a new one.
:::

To use the key, save it with `cytetype login` or `LoginCyteTypeR()`, or pass it as `auth_token`. [Remote servers, SSH and headless machines](../getting-started/sign-in.md#sign-in-remote) walks through it.

In the organization workspace, only admins can create keys. A member who clicks **Create API key** sees a dialog listing the organization's admins to contact.

### The key list

| Column | Meaning |
|---|---|
| Owner | The member who owns the key (shown to organization admins) |
| Name | The key name |
| Key ID | The key's prefix and last four characters. It identifies the key but cannot be used as a key. |
| Status | Active or Disabled |
| Jobs, Clusters | Running and completed jobs submitted with the key, and their clusters, including deleted jobs |
| Quota | The key's cluster limit, or Unlimited |
| Created | The creation date |

The **Active** filter shows active keys only; **All** also shows disabled keys. The **Sort** menu orders keys by date created, name, clusters used, jobs run or quota.

### Rename a key or change its quota

Click the edit icon next to the name or the quota, change the values and click **Save Changes**. Leave the quota empty for no per-key limit. In your personal workspace you can edit your own keys. In the organization workspace, only admins can edit keys, and they can edit any organization key.

### Disable a key

Click **Disable** and confirm. The key stops working at once: uploads and job submissions that use it are rejected with the error code `TOKEN_INACTIVE` (an `AuthenticationError` in Python), and it can no longer fetch the results of private or organization reports. Jobs that are already running still finish; fetch their results with another valid key.

Disabling is permanent. Keys cannot be re-enabled, and they cannot be deleted: disabled keys stay in the list under **All** as a record. Members can disable their own organization keys; admins can disable any organization key.

### The "CyteType CLI" key

Browser sign-in (`cytetype setup` or `SetupCyteTypeR()`) does not always create a new key. It hands out your active personal key named "CyteType CLI", including its secret, and creates that key only when you have no active one. So:

- every computer you sign in from this way shares one key, and each browser sign-in returns its secret again;
- disabling the "CyteType CLI" key here signs all of those computers out, and the next browser sign-in creates a new one ([Switch accounts or servers](../getting-started/sign-in.md#sign-in-switch) explains how to sign in again after a key was rejected).

For a computer or pipeline that should have a key of its own, create one on this page instead ([Remote servers, SSH and headless machines](../getting-started/sign-in.md#sign-in-remote)).

`cytetype logout` and `LogoutCyteTypeR()` only delete the saved copy on that computer. To stop a key from working anywhere, disable it here.

## Account

Open the Account page from your email address at the bottom of the sidebar, or at `https://cytetype.nygen.io/account`. It shows the selected workspace:

- **Profile**: the email address of your account.
- **Workspace**: the workspace name, your role, the license type and, for an organization, the license end date ("No end date" for licenses without one).
- **Usage and limits**:
  - Personal workspace: the jobs remaining in the last 24 hours, 7 days and 30 days, with how many of each you have used, and the total clusters and datasets you have annotated.
  - Organization workspace: the clusters annotated in the current license period, the cluster limit, the datasets processed and the clusters remaining.
- **Appearance**: System, Light or Dark.
- **Support and resources**: **Contact support** (an email to contact@nygen.io), the documentation on GitHub and, in the personal workspace, a link to book a meeting when one is available.
- **Session**: **Sign out** of this browser.

[Usage limits and quotas](limits.md) explains what counts toward each number.

## Organization roles

Organization members are either admins or members. Admin rights cover API keys; they give no extra access to reports.

| Action | Member | Admin |
|---|---|---|
| Manage their own reports (rename, visibility, share, transfer, archive, delete) | Yes | Yes |
| Open other members' Organization and Public reports | Yes | Yes |
| Open other members' private reports | Only if shared with them | Only if shared with them |
| Create organization API keys | No | Yes, for themselves or any active member |
| See organization API keys | Their own | All, with owner and usage |
| Rename keys and change quotas | No | Any organization key |
| Disable keys | Their own | Any organization key |

Admins cannot see members' personal workspaces (personal reports and personal keys) or their Chat conversations, which are private to each user.

The dashboard has no controls for organization membership or roles. A new account whose email domain belongs to an organization with an active license joins it as a member at its first sign-in. To change admins or license settings, contact contact@nygen.io.

## See also

- [Visibility and sharing](sharing.md): what each visibility level and an email share allow.
- [Usage limits and quotas](limits.md): what counts toward the numbers on the Account page.
- [Sign in and API keys](../getting-started/sign-in.md): save a key from the API keys page on another computer.
